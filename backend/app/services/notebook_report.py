"""A report freezes one collection of the same recipes used by individual exports."""
from contextvars import ContextVar
from dataclasses import dataclass, field, replace
import ast
import io
import json
import logging

import pandas as pd

from app.services.notebook_archive import NotebookRecipe, dataframe_bytes, portable_filename
from app.services.notebook_preprocessing import report_recipe_sink, report_replays
from app.services.notebook_snapshots import snapshots
from app.services.notebook_source import (
    NotebookSource, extract_notebook_source, function_definition, merge_source, module_source,
)

LOGGER = logging.getLogger(__name__)
_active = ContextVar('report_notebook', default=None)


@dataclass
class ReportRecipes:
    recipes: list[NotebookRecipe] = field(default_factory=list)
    # A report notebook is offered only when every displayed section was captured.
    omitted: list[str] = field(default_factory=list)

    def append(self, recipe):
        self.recipes.append(recipe)


_SECTION_NAMES = {
    'qc_coverage_plot': 'Coverage', 'qc_boxplot_plot': 'Intensity Boxplot', 'qc_cv_plot': 'Coefficient of Variation',
    'qc_intensity_histogram_plot': 'Intensity Distribution', 'qc_correlation_plot': 'Sample Correlation',
    'qc_pca_plot': 'PCA', 'qc_pca_interactive_html': 'Interactive PCA', 'qc_abundance_plot': 'Abundance',
    'qc_abundance_interactive_html': 'Interactive Abundance', 'qc_peptide_coverage_plot': 'Coverage',
    'peptide_rt_plot': 'Retention Time', 'peptide_modification_plot': 'Modifications',
    'peptide_missed_cleavage_plot': 'Missed Cleavages', 'phosphosite_plot_png': 'Phosphosite Counts',
    'phospho_coverage_png': 'Phosphosite Coverage', 'phospho_sty_png': 'STY Plot', 'completeness_missing_value_plot': 'Missing Value Plot',
    'completeness_missing_value_plot_peptide': 'Missing Value Plot', 'volcano_png': 'Volcano Plot',
    'volcano_html': 'Interactive Volcano Plot', 'qc_coverage_table': 'Coverage summary table',
    'qc_boxplot_table': 'Boxplot summary table', 'qc_cv_table': 'CV summary table',
    'pearson_correlation_png': 'Pearson Correlation', 'pearson_correlation_html': 'Interactive Pearson Correlation',
    'venn_png': 'Venn Diagram', 'single_protein_boxplot_plot': 'Single Protein Boxplot',
    'single_protein_lineplot_plot': 'Single Protein Line Plot', 'anova_png': 'One-way ANOVA',
    'ora_plot_png': 'GO Enrichment (ORA)',
}


def _section_label(renderer, args, options):
    """How the report names this part, for the list of sections the notebook leaves out."""
    name = getattr(renderer, '__name__', str(renderer))
    first = args[0] if args else None
    kind = options.get('kind') or getattr(first, 'kind', None) or (first if isinstance(first, str) else '')
    if name.startswith('peptide_') or name.endswith('_peptide') or name == 'qc_peptide_coverage_plot':
        kind = 'peptide'
    label = _SECTION_NAMES.get(name, name.replace('_', ' '))
    if getattr(first, 'condition1', None):
        label += f' ({first.condition1} vs {first.condition2})'
    return f'{str(kind).capitalize()} — {label}' if kind else label


def _volcano_counts(table):
    """The counts the report's volcano paragraph states, from the volcano table itself."""
    counts = {'upregulated': int((table['significance'] == 'Upregulated').sum()),
              'downregulated': int((table['significance'] == 'Downregulated').sum()),
              'notSignificant': int((table['significance'] == 'Not significant').sum()),
              'unavailable': int((table['significance'] == 'Unavailable').sum()),
              'tested': int(len(table))}
    significant = table[table['significance'].isin(['Upregulated', 'Downregulated'])]
    if len(significant) and 'imputed_values' in table:
        counts['significantWithImputedValues'] = int((significant['imputed_values'].fillna(0) > 0).sum())
    return counts


_VOLCANO_COUNTS_CODE = (
    "\ncounts = functions.functions_report_counts(table)\n"
    "assert counts == json.loads((RESOURCES / 'counts.json').read_text(encoding='utf-8'))\n"
    "Path('reproduced_counts.json').write_text(json.dumps(counts, indent=2), encoding='utf-8', newline='')\n"
    "print(f\"{counts['upregulated']} upregulated, {counts['downregulated']} downregulated, \"\n"
    "      f\"{counts['notSignificant']} not significant, {counts['unavailable']} without a p-value, \"\n"
    "      f\"out of {counts['tested']} tested.\")\n")


def _render(name, *args, **options):
    from app.services import notebook_qc as qc, notebook_pca as pca
    from app.services import notebook_peptide as peptide, notebook_phospho as phospho
    from app.services import notebook_completeness as completeness, notebook_interactive as interactive
    from app.services.notebook_volcano import render_volcano
    from app.services import notebook_comparison as comparison, notebook_single as single
    from app.services.notebook_anova import render_anova
    from app.services.notebook_ora import render_ora

    qc_names = {renderer: plot for plot, renderer in qc.QC_RENDERERS.items()}
    if name in qc_names:
        return qc.render_qc(*args, plot=qc_names[name], **options)
    adapters = {
        'qc_pca_plot': pca.render_pca,
        'qc_pca_interactive_html': lambda **kw: pca.render_pca(interactive=True, **kw),
        'qc_abundance_interactive_html': interactive.render_abundance,
        'qc_peptide_coverage_plot': peptide.render_coverage,
        'peptide_rt_plot': peptide.render_rt,
        'peptide_modification_plot': peptide.render_modification,
        'peptide_missed_cleavage_plot': peptide.render_missed_cleavage,
        'phosphosite_plot_png': phospho.render_counts,
        'phospho_coverage_png': phospho.render_coverage,
        'phospho_sty_png': phospho.render_sty,
        'completeness_missing_value_plot': lambda **kw: completeness.render_completeness(plot='missing-value', **kw),
        'completeness_missing_value_plot_peptide': lambda **kw: completeness.render_peptide_missing_values('peptide', **kw),
        'volcano_png': render_volcano,
        'volcano_html': lambda *a, **kw: render_volcano(*a, interactive=True, **kw),
        'pearson_correlation_png': lambda **kw: comparison.render_comparison('pearson', **kw),
        'pearson_correlation_html': lambda **kw: comparison.render_comparison('pearson', interactive=True, **kw),
        'venn_png': lambda **kw: comparison.render_comparison('venn', **kw),
        'single_protein_boxplot_plot': single.render_box,
        'single_protein_lineplot_plot': single.render_line,
        'anova_png': render_anova,
        'ora_plot_png': render_ora,
    }
    return adapters[name](*args, **options)


def report_plot(renderer, *args, **options):
    collector = _active.get()
    if collector is None:
        return renderer(*args, **options)
    before = len(collector.recipes)
    try:
        content, token = _render(renderer.__name__, *args, **options)
        if not token or len(collector.recipes) != before + 1:
            del collector.recipes[before:]
            collector.omitted.append(_section_label(renderer, args, options))
        else:
            recipe = collector.recipes[-1]
            kind = recipe.provenance.get('kind', '')
            if recipe.provenance.get('plot') == 'volcano':
                table = pd.read_parquet(io.BytesIO(recipe.resources['table.parquet']))
                recipe = replace(recipe,
                    resources={**recipe.resources, 'counts.json': json.dumps(_volcano_counts(table)).encode()},
                    source=extract_notebook_source(_merge_counts(recipe.source.code),
                        (*recipe.source.functions, 'functions_report_counts')),
                    calculation=recipe.calculation + _VOLCANO_COUNTS_CODE)
            collector.recipes[-1] = replace(recipe, plot_name=f'{kind.capitalize()} — {recipe.plot_name}')
        return content
    except Exception:
        # An unavailable export must not take a working HTML report away.
        LOGGER.exception('Report notebook adapter failed: %s', renderer.__name__)
        del collector.recipes[before:]
        # A plot that cannot be drawn either is left out of the report, so the
        # notebook is not missing it.
        content = renderer(*args, **options)
        collector.omitted.append(_section_label(renderer, args, options))
        return content


def _merge_counts(code):
    counts = _definition_renamed(module_source('notebook_report.py'), '_volcano_counts', 'functions_report_counts')
    return merge_source(code, counts)


def _definition_renamed(source, name, new_name):
    return function_definition(source, name).replace(f'def {name}(', f'def {new_name}(', 1)


def report_table(renderer, *args, **options):
    frame = renderer(*args, **options)
    collector = _active.get()
    if collector is None or frame.empty:
        return frame
    plot = {'qc_coverage_table': 'coverage', 'qc_boxplot_table': 'boxplot', 'qc_cv_table': 'cv'}[renderer.__name__]
    kind = options.get('kind', args[0] if args else None)
    for recipe in reversed(collector.recipes):
        if recipe.provenance.get('kind') == kind and recipe.provenance.get('plot') == plot:
            table_bytes = recipe.resources.get('table.parquet')
            if table_bytes is not None and pd.read_parquet(io.BytesIO(table_bytes)).equals(frame):
                return frame  # Already calculated alongside the matching figure.
    before = len(collector.recipes)
    try:
        from app.services.notebook_qc import render_qc
        _, token = render_qc(*args, plot=plot, **options)
        if not token or len(collector.recipes) != before + 1:
            raise ValueError('Missing report table recipe.')
        recipe = collector.recipes[-1]
        pd.testing.assert_frame_equal(pd.read_parquet(io.BytesIO(recipe.resources['table.parquet'])), frame)
        code = recipe.calculation[recipe.calculation.index('table_options ='):]
        roots = tuple(dict.fromkeys(node.attr for node in ast.walk(ast.parse(recipe.preparation + '\n' + code))
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == 'functions'))
        if recipe.replay:
            # The input cell is written at download; it reads and replays the original file.
            roots = tuple(dict.fromkeys((*roots, 'read_original_table', 'replay_id_translations',
                                         'place_imputed_values', recipe.replay['function'])))
        collector.recipes[-1] = replace(recipe, plot_name=f'{kind.capitalize()} {plot} summary table',
            source=extract_notebook_source(recipe.source.code, roots),
            resources={key: value for key, value in recipe.resources.items() if not key.startswith('original.')},
            provenance={**recipe.provenance, 'outputs': ['table']},
            calculation='from IPython.display import display\n' + code)
    except Exception:
        LOGGER.exception('Report summary table notebook unavailable')
        del collector.recipes[before:]
        collector.omitted.append(_section_label(renderer, args, options))
    return frame


def record_sample_mapping(title, mapping):
    collector = _active.get()
    if collector is None:
        return
    try:
        metadata = mapping.drop(columns=['label'])
        source = extract_notebook_source(module_source('functions.py'), ('sample_label_frame',))
        collector.append(NotebookRecipe(plot_name=f'Sample mapping — {title}', source=source,
            settings={}, provenance={'inputReplay': 'captured-metadata', 'outputs': ['table']},
            original_filename='', resources={'input.parquet': dataframe_bytes(metadata),
                'metadata.parquet': dataframe_bytes(metadata), 'table.parquet': dataframe_bytes(mapping)},
            input_explanation='Sample labels are regenerated from the report’s captured annotation metadata.',
            calculation="from IPython.display import display\n"
                "table = functions.sample_label_frame(metadata)\n"
                "if 'custom_label' in metadata:\n"
                "    table['custom_label'] = table['sample'].map(dict(zip(metadata['sample'], metadata['custom_label'])))\n"
                "pd.testing.assert_frame_equal(table, read_resource_table(RESOURCES / 'table.parquet'), check_exact=True)\n"
                "table.to_csv('reproduced_table.csv', index=False)\ndisplay(table)"))
    except Exception:
        LOGGER.exception('Report sample mapping notebook unavailable')
        collector.omitted.append(f'Sample mapping — {title}')


def load_functions(path, name):
    """Separate modules keep section-specific scientific helpers from colliding."""
    import importlib.util
    import sys
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def combine_report(title, recipes, ship_matrices=False, ship_imputed=False):
    """One notebook for the whole report, packaged for the two export options.

    Each section keeps its own code and settings. A data matrix several
    sections share is stored once, an original file is read once, and a
    preprocessing is run once.
    """
    from hashlib import sha256
    from app.services.notebook_archive import (
        BASIS_MATRIX, INPUT_MATRIX, export_explanation, imputed_values_note, protein_side_code, replay_preparation, shipped_matrices,
        side_resource,
    )

    loader = extract_notebook_source(module_source('notebook_report.py'), ('load_functions',))
    imports = tuple(sorted(set(loader.imports).union(*(recipe.source.imports for recipe in recipes))))
    resources = {'input.parquet': dataframe_bytes(pd.DataFrame()),
                 'metadata.parquet': dataframe_bytes(pd.DataFrame())}
    sections, provenance, inputs, recorded, required = [], [], {}, {}, set()
    for index, recipe in enumerate(recipes, 1):
        section = f'section_{index:03d}'
        replay = recipe.replay
        shipped = shipped_matrices(replay, ship_matrices, ship_imputed) if replay else set()
        protein = replay.get('protein') if replay else None
        protein_shipped = shipped_matrices(protein, ship_matrices, ship_imputed) if protein else set()
        matrices, protein_matrices = {}, {}
        for key, value in recipe.resources.items():
            for base, wanted, found in ((INPUT_MATRIX, shipped, matrices), (BASIS_MATRIX, shipped, matrices),
                                        (side_resource('protein', INPUT_MATRIX), protein_shipped, protein_matrices),
                                        (side_resource('protein', BASIS_MATRIX), protein_shipped, protein_matrices)):
                if key == base and (replay if found is matrices else protein):
                    if key.removeprefix('protein_') in wanted:
                        name = f'shared/{sha256(value).hexdigest()[:16]}.parquet'
                        resources[name] = value
                        found[key.removeprefix('protein_')] = f"ROOT_RESOURCES / {json.dumps(name)}"
                    break
            else:
                resources[f'{section}/{key}'] = value
        resources[f'{section}/functions.py'] = recipe.source.code.encode()
        resources[f'{section}/settings.json'] = json.dumps(recipe.settings, allow_nan=False).encode()
        filename = portable_filename(recipe.original_filename)
        # Warnings belong to the screen, never to the report or its export.
        section_provenance = {key: value for key, value in recipe.provenance.items() if key != 'warnings'}
        item = {'section': section, 'name': recipe.plot_name, 'originalFilename': filename,
                'usesImputedInput': recipe.uses_imputed_input, 'provenance': section_provenance}
        provenance.append(item)
        code = (f"RESOURCES = ROOT_RESOURCES / '{section}'\n"
            f"functions = report_functions.load_functions(RESOURCES / 'functions.py', '{section}')\n"
            "settings = json.loads((RESOURCES / 'settings.json').read_text(encoding='utf-8'))\n"
            "metadata = read_resource_table(RESOURCES / 'metadata.parquet')\n"
            f"Path('outputs/{section}').mkdir(parents=True, exist_ok=True)\n")
        preparation, explanation, calculation = recipe.preparation, recipe.input_explanation, recipe.calculation
        if recipe.provenance.get('inputReplay') == 'original-file':
            # A collapsed table is rebuilt from the phosphosite file, the same one
            # the phosphosite sections read.
            key = 'phospho' if replay and replay.get('collapse') else recipe.provenance['kind']
            if key in inputs and inputs[key] != filename:
                raise ValueError('Conflicting original uploads in the report.')
            inputs[key] = filename
            if replay:
                recorded[section] = replay['recorded']
                if not shipped:
                    required.add(key)
                preparation = replay_preparation(replay, shipped, {'section': section, 'path': key, 'matrices': matrices})
                explanation = export_explanation(replay, shipped)
                item['provenance'] = {**section_provenance, 'shippedMatrices': sorted(shipped)}
                if protein:
                    protein_file = portable_filename(protein['recorded']['filename'])
                    if inputs.setdefault('protein', protein_file) != protein_file:
                        raise ValueError('Conflicting original uploads in the report.')
                    recorded[section] = {**replay['recorded'], 'protein': protein['recorded']}
                    if not protein_shipped:
                        required.add('protein')
                    calculation, _ = protein_side_code(calculation, replay, ship_matrices, ship_imputed,
                        {'section': section, 'path': 'protein', 'matrices': protein_matrices})
            else:
                required.add(key)
                code += f'data_path = data_paths[{json.dumps(key)}]\n'
        if recipe.uses_imputed_input:
            if not recipe.imputation_recipe.strip():
                raise ValueError('Imputation history is required for every imputed report section.')
            code += '# The imputed values are the ones CoPYlot drew; no new values are drawn.\n'
            note = imputed_values_note(recipe.imputation_recipe, replay)
            code += ''.join('# ' + line + '\n' for line in note.splitlines())
            item['imputationRecipe'] = note
        code += preparation + '\n' + calculation
        # Existing recipes use fixed output basenames; isolate those, not their calculations.
        code = code.replace("'reproduced", f"'outputs/{section}/reproduced")
        code = code.replace("'recalculated", f"'outputs/{section}/recalculated")
        code = code.replace('"reproduced', f'"outputs/{section}/reproduced')
        sections.append((f'## {index}. {recipe.plot_name}\n\n{explanation}', code))
    preparation = ('report_functions = functions\nROOT_RESOURCES = RESOURCES\n\n'
                   '# The original uploaded files (filenames below). A path marked required must be set;\n'
                   '# an optional one may stay empty, the shipped matrices are then used.\n'
                   'data_paths = {\n')
    for kind, filename in inputs.items():
        value, note = (('Path("YOUR_ORIGINAL_DATA_FILE")', 'required') if kind in required else ('""', 'optional'))
        preparation += f'    {json.dumps(kind)}: {value},  # {note}: {json.dumps(filename)}\n'
    preparation += ('}\n\n'
                    '# A file is read once and a preprocessing is run once for all sections that use it.\n'
                    '_shared = {}\n'
                    'def shared(key, compute):\n'
                    '    if key not in _shared:\n'
                    '        _shared[key] = compute()\n'
                    '    return _shared[key]\n')
    return NotebookRecipe(plot_name=title + ' — Report',
        source=NotebookSource(loader.code, loader.functions, imports), settings={},
        provenance={'inputReplay': 'report-sections', 'sections': provenance, 'originalFiles': inputs,
                    'requiredOriginalFiles': sorted(required), 'recordedInputs': recorded, 'omittedSections': []},
        original_filename='', resources=resources, preparation=preparation, calculation='',
        uses_imputed_input=any(recipe.uses_imputed_input for recipe in recipes),
        imputation_recipe='Sections marked usesImputedInput use the imputed values CoPYlot drew. Exact parameters '
            'are documented in those sections and in the manifest.',
        calculation_sections=tuple(sections), input_explanation=
            'This notebook reproduces the plots, summary tables and volcano counts of the generated report. '
            'The report\u2019s prose (methods paragraphs and any text you wrote) is not regenerated. '
            'Set the original-file paths below, then run from top to bottom. Each section contains its '
            'own calculation code and settings; outputs are saved in separate numbered folders. '
            'An original file is read once and each preprocessing is run once for all sections that use it. '
            'The manifest identifies imputation separately for every section.')


def report_with_notebook(builder, payload):
    from app.services.runtime_cache import runtime_cache_version
    from app.services.session_context import get_session_id
    revision = runtime_cache_version()
    collector = ReportRecipes()
    active_token = _active.set(collector)
    sink_token = report_recipe_sink.set(collector.append)
    replay_token = report_replays.set({})
    try:
        result = builder(payload)
    finally:
        report_replays.reset(replay_token)
        report_recipe_sink.reset(sink_token)
        _active.reset(active_token)
    try:
        from app.services.viewer_policy import notebook_exports_blocked
        if notebook_exports_blocked():
            return result  # Viewers download the report, not its notebook.
        if runtime_cache_version() != revision:
            raise ValueError('The data changed while the report was generated.')
        if collector.omitted:
            raise ValueError('Report sections could not be captured: ' + '; '.join(collector.omitted))
        if collector.recipes:
            title, recipes = payload.title.strip() or 'CoPYlot', tuple(collector.recipes)
            # Packaged again at download for the user's export options.
            recipe = replace(combine_report(title, recipes),
                assemble=lambda matrices, imputed: combine_report(title, recipes, matrices, imputed))
            result.notebookToken = snapshots.capture(get_session_id(), recipe)
    except Exception:
        LOGGER.exception('Full report notebook unavailable')
        result.notebookError = ('The complete report notebook could not be prepared. '
            'No partial notebook was exported. Regenerate the report; if this persists, '
            'check the server log for the affected section.')
    return result
