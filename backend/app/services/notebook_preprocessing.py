"""Verified original-file preprocessing for table-backed notebook exports."""
from dataclasses import replace
from contextvars import ContextVar
from functools import lru_cache, wraps
import hashlib
import io
import inspect
import json
import ast
from types import ModuleType

import pandas as pd

from app.services.notebook_archive import (
    BASIS_MATRIX, IMPUTED_VALUES, NOTEBOOK_HELPERS, dataframe_bytes, place_imputed_values, table_fingerprint,
)
from app.services.notebook_provenance import NotebookProvenanceError
from app.services.notebook_snapshots import PendingRecipe, RetainedReplay, snapshots
from app.services.notebook_source import (
    extract_notebook_source, function_definition, merge_source, module_source, replace_once,
)
from app.services.dataset_reader import adapt_report, rename_cols, sanitize_dataframe
from app.services.annotation_processor import filter_then_normalize
from app.services.contaminant_filter import contaminant_reference
from app.services.id_translation import replay_id_translations
from app.services.functions import filter_localization, inverse_log2_transform_data, raw_zero_to_nan
from app.services.phosprot_tools import transform_phosprot


report_recipe_sink = ContextVar('report_recipe_sink', default=None)


_render_recipe = ContextVar("render_recipe", default=None)
report_replays = ContextVar("report_replays", default=None)


def _render_history_key():
    """Capture editable scientific history, including upstream records used by combined plots."""
    from app.services.annotation_store import get_annotation
    from app.services.dataset_store import get_current_dataset

    history = []
    for kind in ("protein", "phospho", "phosprot"):
        annotation = get_annotation(kind)
        dataset = get_current_dataset(kind)
        history.append((kind, annotation.imputation if annotation else None,
                        annotation.phosprot_lineage if annotation else None,
                        annotation.phosprot_aggregation if annotation else None,
                        dataset.phosprot_lineage if dataset else None,
                        dataset.phosprot_aggregation if dataset else None))
    return history


def cached_notebook_render(render):
    """Reuse the figure and its immutable inputs, never a previous request's token."""
    from app.services.runtime_cache import cached_call, runtime_cache_version
    from app.services.plot_export import plot_notices, notebook_unavailable, note_plot_notice, note_notebook_unavailable
    from app.services.session_context import get_session_id
    signature = inspect.signature(render)

    @wraps(render)
    def wrapped(*args, **kwargs):
        revision = runtime_cache_version()
        capture = notebook_capture_requested()
        arguments = signature.bind(*args, **kwargs)
        arguments.apply_defaults()

        def draw():
            recipes = []
            notices = []
            unavailable = {}
            recipe_token = _render_recipe.set(recipes)
            notices_token = plot_notices.set(notices)
            error_token = notebook_unavailable.set(unavailable)
            try:
                content, _ = render(*args, **kwargs)
            finally:
                _render_recipe.reset(recipe_token)
                plot_notices.reset(notices_token)
                notebook_unavailable.reset(error_token)
            return content, next(iter(recipes), None), notices, unavailable.get("reason")

        key = (capture, arguments.arguments, _render_history_key())
        content, recipe, notices, reason = cached_call(
            render.__module__ + "." + render.__name__, key, draw, copy_result=False)
        for notice in notices:
            note_plot_notice(notice)
        if reason:
            note_notebook_unavailable(reason)
        token = ""
        if capture and recipe is not None and runtime_cache_version() == revision:
            try:
                token = capture_notebook(get_session_id(), recipe)
            except ValueError:
                pass  # capture_notebook already recorded the unavailable reason.
        return content, token
    return wrapped


def notebook_capture_requested():
    from app.services.plot_export import notebook_capture_allowed
    from app.services.viewer_policy import notebook_exports_blocked

    return report_recipe_sink.get() is not None or (notebook_capture_allowed.get() and not notebook_exports_blocked())


def capture_notebook(session, recipe):
    """Snapshot a rendered output for export; '' when this render gets no notebook.

    A plot render only binds the replay (see `bind_original_preprocessing`); it
    is prepared and checked when the ZIP is downloaded. A report prepares every
    section at once, because the report is one deliberate export.
    """
    from copy import deepcopy
    from app.services.runtime_cache import runtime_cache_version
    sink = report_recipe_sink.get()
    if not notebook_capture_requested():
        return ''
    deferred = _render_recipe.get()
    if deferred is not None:
        deferred.append(recipe)
        return "deferred"
    from app.services.plot_export import note_notebook_unavailable
    try:
        revision = runtime_cache_version()
        recipe = deepcopy(recipe)
        finish = bind_original_preprocessing(recipe)
        if sink is not None:
            recipe = finish()
        if runtime_cache_version() != revision:
            raise NotebookProvenanceError('The dataset changed during notebook preparation; render again.')
        if sink is not None:
            sink(recipe)
            return 'report-section'
        return snapshots.capture(session, PendingRecipe(recipe, finish))
    except ValueError as exc:
        # The adapters keep the plot and drop the token; the page says why.
        note_notebook_unavailable(exc)
        raise


def _drawn_values(imputed, basis):
    """The cells imputation filled, as (row label, column, value); None when the
    matrix is anything other than its basis with missing cells filled."""
    import numpy as np
    if basis is None or imputed is None or not imputed.index.equals(basis.index) or list(imputed.columns) != list(basis.columns):
        return None
    rows, columns, values = [], [], []
    for position, column in enumerate(imputed.columns):
        after, before = imputed.iloc[:, position], basis.iloc[:, position]
        if pd.api.types.is_numeric_dtype(after) and pd.api.types.is_numeric_dtype(before):
            a, b = after.to_numpy(dtype=float), before.to_numpy(dtype=float)
            differs = ~((a == b) | (np.isnan(a) & np.isnan(b)))
            if (differs & ~np.isnan(b)).any() or not np.isfinite(a[differs]).all():
                return None  # An observed value changed: not a fill of missing cells.
            found = np.flatnonzero(differs)
            rows.extend(imputed.index[found].tolist())
            columns.extend([str(column)] * len(found))
            values.extend(a[found].tolist())
        elif not after.astype(object).where(after.notna(), None).equals(before.astype(object).where(before.notna(), None)):
            return None
    return pd.DataFrame({'row': rows, 'column': columns, 'value': values})


def read_original_table(path, select_column=None, peptide=False):
    """Use the upload reader's report adapters, column aliases and missing identifier handling."""
    from pathlib import Path
    path = Path(path)
    extension = path.suffix.lower()
    if extension in {'.csv', '.tsv', '.txt'}:
        chunks = list(pd.read_csv(path, sep=',' if extension == '.csv' else '\t', chunksize=50_000, usecols=select_column))
        frame = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    elif extension == '.xlsx':
        import openpyxl  # noqa: F401 -- pandas reads the original workbook with it.
        frame = pd.read_excel(path, usecols=select_column)
    elif extension == '.parquet':
        import pyarrow.parquet as pq
        parquet = pq.ParquetFile(path)
        columns = [c for c in parquet.schema_arrow.names if select_column(c)] if select_column else None
        batches = [batch.to_pandas() for batch in parquet.iter_batches(batch_size=50_000, columns=columns)]
        frame = (pd.concat(batches, ignore_index=True) if len(batches) > 1 else batches[0]) if batches else parquet.read(columns=columns).to_pandas()
    else:
        raise ValueError('Use the original CSV, TSV, TXT, XLSX or Parquet upload.')
    return sanitize_dataframe(rename_cols(adapt_report(frame, peptide)))


def replay_preprocessing(raw, metadata, config):
    raw = replay_id_translations(raw, config.get('id_translations', []))
    samples = metadata['sample'].tolist()
    if config.get('vsn_runtime') and vsn_runtime_version() != config['vsn_runtime']:
        raise ValueError('Use the recorded R and Bioconductor vsn versions for this notebook.')
    if config['source'] == 'original_source_rows':
        data = raw.copy()
    elif config['source'] == 'source_data':
        data = raw.copy() if config['is_log2'] else raw_zero_to_nan(raw, samples)
    else:
        # The app's own order: the filter decides the rows, the normalization is
        # fitted on them and applied to every row (filter_then_normalize).
        filtering = config['filter']
        log2_data, filtered_data = filter_then_normalize(raw, metadata, config['is_log2'], config['normalization'],
            min_present=filtering['minPresent'], filter_mode=filtering['mode'],
            filter_method=filtering['method'], min_percentage=filtering['minPercentage'],
            normalization_fit=filtering.get('normalizationFit', 'filtered_rows'),
            contaminants=filtering.get('contaminants'), contaminant_snapshot=config.get('contaminant_reference'))
        data = filtered_data if config['source'] == 'filtered' else log2_data
    if config['kind'] == 'phospho' and config.get('localization_cutoff') is not None:
        data = filter_localization(data, config['localization_cutoff'])
    if config['linear_output']:
        data = inverse_log2_transform_data(data, samples).reindex(columns=data.columns)
    return data


def replay_collapsed_input(raw, metadata, config, upstream_imputed_values=None):
    upstream_metadata = pd.DataFrame(config['upstream_metadata'])
    observed = sites = replay_preprocessing(raw, upstream_metadata, config['upstream'])
    if upstream_imputed_values is not None:
        # The collapse used imputed phosphosites: the drawn values are placed, never drawn again.
        # Site coverage still counts what was measured.
        sites = place_imputed_values(observed, upstream_imputed_values)
    collapsed = transform_phosprot(sites, upstream_metadata, mode=config['aggregation'],
        data_scale='log2', coverage_data=observed)
    return replay_preprocessing(collapsed, metadata, config)


def vsn_runtime_version():
    import subprocess
    from app.services.r_runtime import detect_r_runtime
    executable = detect_r_runtime().get('executable')
    if not executable:
        raise ValueError('Install R and Bioconductor vsn to replay this normalization.')
    try:
        result = subprocess.run([str(executable), '-e',
            'cat(paste(getRversion(), as.character(packageVersion("vsn")), sep="|"))'],
            capture_output=True, text=True, check=True, timeout=60,
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError('R and Bioconductor vsn must be available to replay this normalization.') from exc
    return result.stdout.strip()


def limma_runtime_version():
    """'Rversion|limmaversion', recorded with a moderated export and checked on replay."""
    import subprocess
    from app.services.r_runtime import detect_r_runtime
    executable = detect_r_runtime().get('executable')
    if not executable:
        raise ValueError('Install R and Bioconductor limma to replay this moderated test.')
    try:
        result = subprocess.run([str(executable), '-e',
            'cat(paste(getRversion(), as.character(packageVersion("limma")), sep="|"))'],
            capture_output=True, text=True, check=True, timeout=60,
            creationflags=getattr(subprocess, 'CREATE_NO_WINDOW', 0))
    except (OSError, subprocess.SubprocessError) as exc:
        raise ValueError('R and Bioconductor limma must be available to replay this moderated test.') from exc
    return result.stdout.strip()


@lru_cache(maxsize=4)
def preprocessing_source(collapse=False, vsn=False):
    functions = module_source('functions.py')
    functions = replace_once(functions, '        from app.services.r_vsn import vsn_normalize_linear', '')
    replay = function_definition(module_source('notebook_preprocessing.py'), 'replay_preprocessing')
    if not vsn:
        start = functions.index('    elif resolved == "vsn":', functions.index('def normalize_data('))
        end = functions.index('\n    else:\n        raise ValueError("Normalization must', start)
        functions = functions[:start] + functions[end:]
        replay = replace_once(replay,
            "    if config.get('vsn_runtime') and vsn_runtime_version() != config['vsn_runtime']:\n"
            "        raise ValueError('Use the recorded R and Bioconductor vsn versions for this notebook.')\n", '')
    code = merge_source(
        functions,
        extract_notebook_source(module_source('dataset_reader.py'), ('adapt_report', 'rename_cols', 'sanitize_dataframe')).code,
        function_definition(module_source('annotation_processor.py'), 'prepare_normalized_log2_data'),
        function_definition(module_source('annotation_processor.py'), 'filter_then_normalize').replace(
            '        reference = contaminant_reference()',
            '        raise ValueError("Missing contaminant snapshot in processing.json.")'),
        extract_notebook_source(module_source('contaminant_filter.py'), ('contaminant_mask',)).code,
        function_definition(module_source('notebook_preprocessing.py'), 'read_original_table'),
        extract_notebook_source(module_source('id_translation.py'), ('replay_id_translations',)).code,
        function_definition(NOTEBOOK_HELPERS, 'place_imputed_values'),
        replay)
    if vsn:
        portable_runtime = ('import shutil\n\ndef detect_r_runtime():\n'
            '    executable = shutil.which("Rscript")\n'
            '    return {"available": bool(executable), "executable": executable}\n')
        vsn_code = module_source('r_vsn.py').replace('from app.services.r_runtime import detect_r_runtime', '')
        version_code = function_definition(module_source('notebook_preprocessing.py'), 'vsn_runtime_version').replace(
            '    from app.services.r_runtime import detect_r_runtime\n', '')
        code = merge_source(code, vsn_code, portable_runtime, version_code)
    root = 'replay_preprocessing'
    if collapse:
        collapse_code = module_source('phosprot_tools.py').replace('PhosprotAggregationMode', 'str')
        # Resolve the original collapse's numeric helpers from the same scientific source.
        nodes = [node for node in ast.parse(collapse_code).body if not isinstance(node, ast.ImportFrom)
                 or not (node.module or '').startswith('app.')]
        collapse_code = '\n\n'.join(ast.get_source_segment(collapse_code, node) for node in nodes)
        code = merge_source(code, collapse_code,
            function_definition(module_source('notebook_preprocessing.py'), 'replay_collapsed_input'))
        root = 'replay_collapsed_input'
    return extract_notebook_source(code, ('read_original_table', 'replay_id_translations', 'place_imputed_values', root))


def _imputation_settings(record):
    """The arguments imputation was run with, or None when the record is incomplete."""
    if record is None or not record.applied:
        return None
    needed = [record.seed, record.sampleWise]
    needed += [record.downShift, record.width] if record.mode == 'numerical' else [record.qValue, record.adjustStd]
    if any(value is None for value in needed):
        return None
    return dict(method=record.mode, seed=record.seed, sample_wise=record.sampleWise,
                q=0.01 if record.qValue is None else record.qValue,
                adj_std=0.3 if record.adjustStd is None else record.adjustStd,
                down_shift=1.8 if record.downShift is None else record.downShift,
                width=0.3 if record.width is None else record.width)


def _processing_config(kind, annotation, dataset, source):
    from app.services.phospho_localization import get_localization_cutoff

    config = {'kind': kind, 'source': source, 'is_log2': annotation.is_log2_transformed,
        'normalization': annotation.filter_config.normalization,
        'filter': annotation.filter_config.model_dump(), 'localization_cutoff': get_localization_cutoff(),
        'linear_output': False, 'id_translations': [step.model_dump() for step in dataset.id_translations]}
    if config['filter'].get('contaminants', {}).get('enabled'):
        config['contaminant_reference'] = contaminant_reference()
    return config


def _bind_protein_side(recipe, source):
    """The protein table of a protein-corrected site comparison, read from live state now.

    It is a second dataset with its own upload and preprocessing. None when it
    has no replay; its matrix then stays frozen in the export, as before.
    """
    from app.services.annotation_store import get_annotation
    from app.services.dataset_store import get_current_dataset

    annotation, dataset = get_annotation('protein'), get_current_dataset('protein')
    if annotation is None or dataset is None:
        return None
    if dataset.id_translations and not annotation.source_data.equals(dataset.frame):
        return None
    imputed = source == 'imputed'
    if imputed:
        if annotation.imputed_data is None or annotation.imputed_data.empty:
            return None
        source = 'filtered' if annotation.filtered_data is not None else 'data'
    if source not in {'data', 'filtered'}:
        return None
    return {'config': _processing_config('protein', annotation, dataset, source), 'raw': annotation.source_data,
            'metadata': annotation.metadata[['sample', 'condition']].copy(), 'imputed': imputed,
            'filename': dataset.filename}


def _finish_protein_side(side, preprocessing, resources):
    """Check the protein replay and add what it needs to `resources`; the replay record or None."""
    config, raw, metadata = side['config'], side['raw'], side['metadata']
    expected = pd.read_parquet(io.BytesIO(resources['protein_input.parquet']))
    try:
        module = ModuleType('_copylot_protein_preprocessing_validation')
        exec(preprocessing.code, module.__dict__)
        replayed = module.replay_preprocessing(raw, metadata, {**config, 'id_translations': []})
        projected = replayed.loc[expected.index, expected.columns]
        drawn = _drawn_values(expected, projected) if side['imputed'] else None
        if side['imputed'] and drawn is None:
            return None
        if not side['imputed']:
            pd.testing.assert_frame_equal(projected, expected, check_exact=True)
    except Exception:
        return None
    resources['protein_processing_metadata.parquet'] = dataframe_bytes(metadata)
    resources['protein_processing.json'] = json.dumps(config, allow_nan=False).encode()
    if not replayed.index.equals(expected.index):
        resources['protein_input_rows.parquet'] = dataframe_bytes(pd.DataFrame(index=expected.index))
    if drawn is not None:
        resources['protein_' + IMPUTED_VALUES] = dataframe_bytes(drawn)
        resources['protein_' + BASIS_MATRIX] = dataframe_bytes(projected)
    return {'function': 'replay_preprocessing', 'collapse': False, 'imputed': drawn is not None, 'callExtra': '',
            'recorded': {
                'replayKey': hashlib.sha256(json.dumps(['replay_preprocessing', config, metadata.to_dict(orient='list')],
                    sort_keys=True, default=str).encode()).hexdigest()[:16],
                'filename': side['filename'], 'sourceFingerprint': table_fingerprint(raw),
                'inputFingerprint': table_fingerprint(expected, index=True),
                'columns': [str(column) for column in expected.columns]}}


def bind_original_preprocessing(recipe):
    """Read everything the replay needs from live state now; return the replay for later.

    Rendering a plot only binds: the processing settings, the metadata and a
    reference to the upload the annotation was built from. The expensive part --
    running the recorded preprocessing again and checking it reproduces the
    plot's input exactly, which for VSN means fitting it in R -- runs when the
    notebook is downloaded, so a plot nobody exports costs nothing extra. Nothing
    here is read from live state at download time.
    """
    from app.services.annotation_store import get_annotation
    from app.services.dataset_store import get_current_dataset
    from app.services.phospho_localization import get_localization_cutoff

    if recipe.provenance.get('inputReplay') in {'original-file', 'submitted-gene-list'}:
        return lambda: recipe
    kind = recipe.provenance.get('kind')
    if kind in {'peptide', 'precursor'}:
        # Deferred on purpose: the peptide binder needs the peptide recipes'
        # sources, which live beside their renderers, and those renderers
        # capture through this module. The edge is real, not a convenience.
        from app.services.notebook_peptide_preprocessing import bind_peptide_preprocessing
        return bind_peptide_preprocessing(recipe)
    if kind not in {'protein', 'phospho', 'phosprot'}:
        raise NotebookProvenanceError('This output has no original-file replay recipe.')
    annotation, dataset = get_annotation(kind), get_current_dataset(kind)
    if annotation is None or dataset is None:
        raise NotebookProvenanceError('Original-file processing history is unavailable.')
    if dataset.id_translations and not annotation.source_data.equals(dataset.frame):
        raise NotebookProvenanceError('Regenerate annotations after identifier translation before exporting a notebook.')
    source = recipe.provenance.get('sourceUsed', recipe.provenance.get('dataSource', 'data'))
    source_requested = source
    imputation_settings = None
    imputed_input = source == 'imputed'
    if imputed_input:
        # The deterministic part of an imputed matrix is the matrix imputation
        # filled; the drawn values are carried and placed into it.
        if annotation.imputed_data is None or annotation.imputed_data.empty:
            return lambda: recipe
        source = 'filtered' if annotation.filtered_data is not None else 'data'
        imputation_settings = _imputation_settings(annotation.imputation)
    elif recipe.uses_imputed_input and kind != 'phosprot':
        return lambda: recipe  # An imputed input this module cannot rebuild stays frozen.
    if source not in {'data', 'filtered', 'source_data', 'original_source_rows'}:
        raise NotebookProvenanceError('This source has no original-file replay recipe yet.')
    config = {'kind': kind, 'source': source, 'is_log2': annotation.is_log2_transformed,
        'normalization': annotation.filter_config.normalization,
        'filter': annotation.filter_config.model_dump(), 'localization_cutoff': get_localization_cutoff(),
        'linear_output': recipe.provenance.get('inputScale') in {'linear', 'linear ratios'},
        'id_translations': [step.model_dump() for step in dataset.id_translations]}
    if recipe.provenance.get('linearQcUnnormalized'):
        # CV and abundance after VSN read the unnormalized data (VSN is not linear).
        config['normalization'] = 'none'
    original_filename = recipe.original_filename
    raw = annotation.source_data
    collapse = kind == 'phosprot' and annotation.phosprot_lineage is not None
    upstream_imputed, upstream_values = False, None
    if collapse:
        lineage = annotation.phosprot_lineage
        upstream_annotation, upstream_dataset = get_annotation('phospho'), get_current_dataset('phospho')
        if upstream_annotation is None or upstream_dataset is None:
            raise NotebookProvenanceError('The original phosphosite upload is needed to replay this collapse.')
        if upstream_dataset.id_translations and not upstream_annotation.source_data.equals(upstream_dataset.frame):
            raise NotebookProvenanceError('Upstream identifier translation is newer than its annotation.')
        config.update(aggregation=lineage.aggregation,
            upstream_metadata=upstream_annotation.metadata[['sample', 'condition']].to_dict(orient='records'),
            upstream={'kind': 'phospho', 'source': 'data' if lineage.source == 'log2' else lineage.source,
                'is_log2': upstream_annotation.is_log2_transformed, 'normalization': lineage.filter.normalization,
                'filter': lineage.filter.model_dump(), 'localization_cutoff': lineage.localizationCutoff,
                'linear_output': False,
                'id_translations': [step.model_dump() for step in upstream_dataset.id_translations]})
        raw = upstream_annotation.source_data
        original_filename = upstream_dataset.filename
        upstream_imputed = lineage.source == 'imputed'
        if upstream_imputed:
            config['upstream']['source'] = 'filtered' if upstream_annotation.filtered_data is not None else 'data'
            upstream_values = _drawn_values(
                upstream_annotation.imputed_data,
                upstream_annotation.filtered_data if upstream_annotation.filtered_data is not None else upstream_annotation.log2_data)
            if upstream_values is None:
                return lambda: recipe
    for processing_config in [config] + ([config['upstream']] if collapse else []):
        if processing_config['filter'].get('contaminants', {}).get('enabled'):
            processing_config['contaminant_reference'] = contaminant_reference()
    processing_metadata = annotation.metadata[['sample', 'condition']].copy()
    protein_side = _bind_protein_side(recipe, source_requested) if 'protein_input.parquet' in recipe.resources else None

    def finish():
        vsn = config['normalization'] == 'vsn' and source in {'data', 'filtered'}
        if collapse:
            vsn = vsn or config['upstream']['normalization'] == 'vsn'
        if protein_side is not None and protein_side['config']['normalization'] == 'vsn':
            vsn = True
        if vsn:
            runtime = vsn_runtime_version()
            config['vsn_runtime'] = runtime
            if collapse and config['upstream']['normalization'] == 'vsn':
                config['upstream']['vsn_runtime'] = runtime
            if protein_side is not None and protein_side['config']['normalization'] == 'vsn':
                protein_side['config']['vsn_runtime'] = runtime
        preprocessing = preprocessing_source(collapse, vsn)
        expected = pd.read_parquet(io.BytesIO(recipe.resources['input.parquet']))
        replay_function = 'replay_collapsed_input' if collapse else 'replay_preprocessing'
        # Annotation inputs already contain ID translations. The notebook's file does not.
        validation_config = {**config, 'id_translations': []} if not collapse else {
            **config, 'upstream': {**config['upstream'], 'id_translations': []}}
        memo = report_replays.get()
        key = (id(raw), replay_function, json.dumps(validation_config, sort_keys=True), id(upstream_values),
               tuple(processing_metadata.columns), tuple(processing_metadata.itertuples(index=False, name=None))) if memo is not None else None
        if memo is not None and key in memo:
            _, replayed = memo[key]
        else:
            module = ModuleType('_copylot_preprocessing_validation')
            exec(preprocessing.code, module.__dict__)
            extra = {'upstream_imputed_values': upstream_values} if upstream_imputed else {}
            replayed = getattr(module, replay_function)(raw, processing_metadata, validation_config, **extra)
            if memo is not None:
                memo[key] = (raw, replayed)  # Retain the source identity until the report is finished.
        drawn = None
        try:
            projected = replayed.loc[expected.index, expected.columns]
            if imputed_input:
                drawn = _drawn_values(expected, projected)
                if drawn is None:
                    raise AssertionError('The imputed matrix is not its basis with the missing cells filled.')
            else:
                pd.testing.assert_frame_equal(projected, expected, check_exact=True)
        except (AssertionError, KeyError) as exc:
            if imputed_input or upstream_imputed:
                # Imputation that no longer belongs to the current matrix cannot be
                # rebuilt from the file. The export keeps the frozen matrix, as before.
                return recipe
            raise NotebookProvenanceError('Recorded preprocessing does not reproduce this plot input; regenerate annotations.') from exc
        show_imputation = drawn is not None and imputation_settings is not None
        merged = extract_notebook_source(
            merge_source(recipe.source.code, preprocessing.code, *(
                [extract_notebook_source(module_source('functions.py'), ('impute_values_with_diagnostics',)).code]
                if show_imputation else [])),
            tuple(dict.fromkeys((*recipe.source.functions, 'read_original_table', 'replay_id_translations',
                                 'place_imputed_values', replay_function,
                                 *(('impute_values_with_diagnostics',) if show_imputation else ())))))
        resources = {**recipe.resources,
            'processing_metadata.parquet': dataframe_bytes(processing_metadata),
            'processing.json': json.dumps(config, allow_nan=False).encode()}
        if not replayed.index.equals(expected.index):
            resources['input_rows.parquet'] = dataframe_bytes(pd.DataFrame(index=expected.index))
        if drawn is not None:
            resources[IMPUTED_VALUES] = dataframe_bytes(drawn)
            resources[BASIS_MATRIX] = dataframe_bytes(projected)
        if upstream_imputed:
            resources['upstream_imputed_values.parquet'] = dataframe_bytes(upstream_values)
        protein_replay = _finish_protein_side(protein_side, preprocessing, resources) if protein_side is not None else None
        runtime_note = ''
        if vsn:
            r_version, vsn_version = runtime.split('|')
            runtime_note = f' Install R {r_version} and Bioconductor vsn {vsn_version}, with Rscript on PATH; these versions are checked before normalization.'
        if recipe.settings.get('limma_runtime'):
            r_version, limma_version = recipe.settings['limma_runtime'].split('|')
            runtime_note += (f' The moderated t-test runs in R {r_version} with Bioconductor limma {limma_version}, '
                'with Rscript on PATH; these versions are checked before the test.')
        projection_note = (' The displayed projection coordinates and cluster labels are retained for exact '
            'orientation; the plotting code reuses them after verifying the rebuilt analysis matrix.'
            if recipe.provenance.get('plot') == 'pca' else '')
        if 'protein_input.parquet' in recipe.resources:
            # The protein side of a protein-corrected site table is a second
            # dataset with its own preprocessing and its own original file.
            projection_note += (" Each site's change is corrected for its protein's change: the protein matrix "
                'of the same package and data source is tested with the same comparison and subtracted with the '
                'standard errors combined. ' + ('It is rebuilt from the protein file like the main input.'
                if protein_replay else 'That matrix is carried as exported, not rebuilt from the protein upload.'))
        if upstream_imputed:
            projection_note += (' The collapse used imputed phosphosites: the values CoPYlot drew are shipped '
                'and placed into the phosphosite matrix before it is collapsed.')
        replay = {
            'function': replay_function, 'collapse': collapse,
            'imputed': True if drawn is not None else ('upstream' if upstream_imputed else False),
            'callExtra': (", upstream_imputed_values=read_resource_table(RESOURCES / 'upstream_imputed_values.parquet')"
                          if upstream_imputed else ''),
            'notes': runtime_note + projection_note,
            'protein': protein_replay,
            # The call that drew the imputed values, shown and run in the notebook.
            'imputation': imputation_settings if show_imputation else None,
            'linear': bool(config['linear_output']),
            'recorded': {
                'replayKey': hashlib.sha256(json.dumps([replay_function, validation_config,
                    processing_metadata.to_dict(orient='list'),
                    table_fingerprint(upstream_values) if upstream_imputed else None],
                    sort_keys=True, default=str).encode()).hexdigest()[:16],
                'filename': original_filename,
                'sourceFingerprint': table_fingerprint(raw),
                'inputFingerprint': table_fingerprint(expected, index=True),
                'columns': [str(column) for column in expected.columns],
            },
        }
        return replace(recipe, source=merged, resources=resources, original_filename=original_filename, replay=replay,
            provenance={**recipe.provenance, 'inputReplay': 'original-file', 'preprocessing': config})

    return RetainedReplay(finish, (raw, processing_metadata) + ((protein_side['raw'],) if protein_side else ()))
