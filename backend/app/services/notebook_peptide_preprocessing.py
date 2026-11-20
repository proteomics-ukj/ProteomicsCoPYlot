"""Original long-format peptide imports for portable plot recipes."""
from dataclasses import replace
from hashlib import sha256
import io
import json
import sys
from types import ModuleType

import pandas as pd

from app.services.functions import raw_zero_to_nan
from app.services.notebook_archive import dataframe_bytes, table_fingerprint
from app.services.notebook_snapshots import RetainedReplay
from app.services.notebook_source import extract_notebook_source, function_definition, merge_source, module_source


def replay_peptide_matrix(raw, config):
    quantity = config['quantity_column']
    if not pd.api.types.is_numeric_dtype(raw[quantity]):
        raw = raw.assign(**{quantity: pd.to_numeric(raw[quantity], errors="coerce")})
    matrix = raw.pivot_table(index=config['feature_column'], columns=config['file_column'],
        values=config['quantity_column'], aggfunc='max')
    rename = {name: config['runs'][sha256(str(name).encode()).hexdigest()]
        for name in matrix.columns if sha256(str(name).encode()).hexdigest() in config['runs']}
    matrix = matrix[list(rename)].rename(columns=rename)[config['columns']]
    matrix.columns.name = config['columns_name']
    return raw_zero_to_nan(matrix, list(matrix.columns))


def bind_peptide_preprocessing(recipe):
    """Read the run mapping from live state now; check the replay when the notebook is downloaded."""
    from app.services import peptide_tools, plot_images, completeness_tools
    from app.services.notebook_peptide import missed_cleavage_source, modification_source
    from app.services.notebook_preprocessing import preprocessing_source
    from app.services.notebook_provenance import NotebookProvenanceError

    plot = recipe.provenance.get('plot', recipe.provenance.get('output'))
    uses_detection_matrix = plot in {'coverage', 'completeness/missing-value'}
    raw = peptide_tools.get_peptide_frame(detection_only=True) if uses_detection_matrix else peptide_tools.get_peptide_frame()
    meta = None
    expected = pd.read_parquet(io.BytesIO(recipe.resources['input.parquet']))
    loader = extract_notebook_source(preprocessing_source().code, ('read_original_table',)).code
    extra = ''
    resources = dict(recipe.resources)
    config = {}
    read, read_key, shipped_lines = "functions.read_original_table(original, peptide=True)", "rows", []
    if uses_detection_matrix:
        pivot, runs = plot_images.peptide_detection_matrix(recipe.provenance['kind'])
        file_column = runs.columns[0]
        if plot == 'coverage':
            meta = plot_images._peptide_coverage_metadata(runs, include_id=recipe.settings['include_id'])
            mapping = peptide_tools._sample_rename_map(pivot.columns,
                meta.rename(columns={'label': 'new_sample'}), recipe.provenance['featureColumn'])
        else:
            runs = [run for run in completeness_tools._peptide_sample_columns(runs) if run in pivot.columns]
            # The existing export disambiguates basenames in this same run order.
            mapping = dict(zip(runs, expected.columns, strict=True))
        config = {'file_column': file_column, 'feature_column': recipe.provenance['featureColumn'],
            'quantity_column': recipe.provenance['quantityColumn'], 'columns': expected.columns.tolist(),
            'columns_name': expected.columns.name,
            'runs': {sha256(str(run).encode()).hexdigest(): label for run, label in mapping.items()}}
        extra = merge_source('from hashlib import sha256',
            function_definition(module_source('notebook_peptide_preprocessing.py'), 'replay_peptide_matrix'),
            extract_notebook_source(module_source('functions.py'), ('raw_zero_to_nan',)).code)
        preparation = "data = functions.replay_peptide_matrix(raw, processing)\n"
        extra = merge_source(extra, extract_notebook_source(module_source('peptide_tools.py'), ('_detection_column',)).code)
        roots = ('replay_peptide_matrix', '_detection_column')
        read, read_key = "functions.read_original_table(original, functions._detection_column, True)", "detection"
    elif plot in {'missed_cleavages', 'modifications'}:
        meta = peptide_tools._prepare_metadata(raw, recipe.settings['include_id'])
        file_column = peptide_tools._file_name_column(raw)
        mapping = peptide_tools._sample_rename_map(raw[file_column].dropna().unique(), meta, file_column)
        config = {'file_column': file_column,
            'runs': {sha256(str(run).encode()).hexdigest(): label for run, label in mapping.items()}}
        meta['sample'] = meta['new_sample']
        resources['processing_metadata.parquet'] = dataframe_bytes(meta)
        full_source = missed_cleavage_source() if plot == 'missed_cleavages' else modification_source()
        extra = extract_notebook_source(full_source.code, ('prepare_input',)).code
        preparation = "from hashlib import sha256\n"
        preparation += "raw = raw.copy()\n"
        preparation += "raw[processing['file_column']] = raw[processing['file_column']].map(\n"
        preparation += "    lambda name: processing['runs'].get(sha256(str(name).encode()).hexdigest(), sha256(str(name).encode()).hexdigest()))\n"
        preparation += "processing_metadata = read_resource_table(RESOURCES / 'processing_metadata.parquet')\n"
        preparation += ("data, sequence_column = functions.prepare_input(raw, processing_metadata)\n"
            if plot == 'missed_cleavages' else "data = functions.prepare_input(raw, processing_metadata)\n")
        roots = ('prepare_input',)
    elif plot in {'retention_time', 'sequence_coverage'}:
        preparation = 'data = raw\n'
        roots = ()
    else:
        raise NotebookProvenanceError('This peptide output has no original-file replay recipe.')
    source = extract_notebook_source(merge_source(recipe.source.code, loader, extra),
        tuple(dict.fromkeys((*recipe.source.functions, 'read_original_table', *roots))))
    resources['processing.json'] = json.dumps(config, allow_nan=False).encode()

    def finish():
        # Check the same generated processing against the actual long-format input.
        module = ModuleType('_copylot_peptide_replay_check')
        sys.modules[module.__name__] = module
        exec(source.code, module.__dict__)
        if uses_detection_matrix:
            replayed = module.replay_peptide_matrix(raw, config)
        elif plot in {'missed_cleavages', 'modifications'}:
            renamed = raw.copy()
            renamed[file_column] = renamed[file_column].map(lambda name: config['runs'].get(
                sha256(str(name).encode()).hexdigest(), sha256(str(name).encode()).hexdigest()))
            replayed = module.prepare_input(renamed, meta)
            if plot == 'missed_cleavages':
                replayed, sequence_column = replayed
        elif roots:
            replayed = module.prepare_points(raw)
        else:
            replayed = raw
        recorded = {'filename': recipe.original_filename, 'replayKey': '',
                    'sourceFingerprint': table_fingerprint(raw),
                    'inputFingerprint': table_fingerprint(expected, index=True),
                    'columns': expected.columns.tolist()}
        if plot == 'missed_cleavages':
            recorded['sequenceColumn'] = sequence_column
            shipped_lines.append("sequence_column = recorded['sequenceColumn']")
        try:
            pd.testing.assert_frame_equal(replayed.loc[expected.index, expected.columns], expected, check_exact=True)
        except (AssertionError, KeyError) as exc:
            raise NotebookProvenanceError('Original peptide processing does not reproduce the displayed input.') from exc
        final = dict(resources)
        if not replayed.index.equals(expected.index):
            final['input_rows.parquet'] = dataframe_bytes(pd.DataFrame(index=expected.index))
        return replace(recipe, source=source, resources=final,
            provenance={**recipe.provenance, 'inputReplay': 'original-file', 'preprocessing': {'kind': 'peptide'}},
            replay={'function': '', 'collapse': False, 'imputed': False, 'recorded': recorded,
                    'custom': {'read': read, 'readKey': read_key, 'body': preparation.splitlines(), 'shipped': shipped_lines},
                    'notes': ' This is the peptide report: the notebook uses its unfiltered rows and repeats loading, '
                             'aggregation and run mapping. Run matching uses hashes, so original filesystem paths '
                             'are not written into the ZIP.'})

    return RetainedReplay(finish, tuple(frame for frame in (raw, expected, meta) if frame is not None))
