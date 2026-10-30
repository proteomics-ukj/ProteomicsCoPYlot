"""Pre-refactor report fingerprints; renderers are deterministic parameter probes."""
import hashlib
import html
import importlib
import json
import re
from datetime import datetime
from types import SimpleNamespace

import pandas as pd
import pytest

from app.schemas.summary import SummaryReportContext, SummaryReportRequest, SummaryVolcanoEntry
from app.services.html_security import GENERATED_HTML_META

report = importlib.import_module('app.services.report_function')
# 2026-09-16: every case that renders a projection moved, because the report now
# asks `qc_pca_feature_summary` how much of the matrix reached the PCA. The
# rendered HTML is byte-identical in all 24 cases -- under this harness the
# summary is one of the stubbed `qc_*` probes, so the completeness sentence
# resolves to nothing; only the recorded call list gained the extra name. The
# sentence itself is covered by test_report_processing_provenance.py.
# 2026-09-17: the six `peptide,...` cases moved, because the summary report now
# passes the page's `Cleavage Rule` through to `peptide_missed_cleavage_plot`
# instead of always rendering under trypsin. Only that one call gained an
# argument (`protease`); under this harness the probe's arguments are what the
# figure's bytes are, so those cases' HTML moves with it. The pass-through
# itself is covered by test_missed_cleavage_reporting.py.
# Layout snapshots use a fixed version; test_report_version_and_ksea_filename
# independently verifies that real reports print the current application version.
# 2026-09-21: the coverage plot's total line became a switch of its own
# (`referenceLine` in the coverage settings), so the report now forwards
# `reference_line=True` to the two coverage renderers. Reverting just that
# argument reproduces the previous hashes; it is the whole difference.
# 2026-09-24: DBSCAN's radius became relative (`dbscan_eps_factor`, x median
# sample distance) instead of an absolute `dbscan_eps`; the PCA renderer
# argument was renamed accordingly. Later the same day the report's projection
# became PCA-only: `method` is always "PCA" and the UMAP/t-SNE arguments are no
# longer forwarded to the renderer. Those renderer arguments are the whole
# difference.
# 2026-09-30: random seeds left the report text (imputation and k-means lines;
# checked by test_report_processing_provenance and test_report_pca_is_always_pca),
# the phospho coverage plots now receive the Phospho-specific page's own
# settings (include-ID, conditions, colours) and state the fixed Class I
# threshold. Those are the whole difference; verified by diffing the rendered
# HTML and renderer calls of every case against the previous code.
# 2026-10-05: the three `failure` cases moved, because warnings no longer go
# into the report: the "Report Warnings" section is gone, the failed figure's
# note no longer points to it, and a warning reads "<figure> could not be
# rendered: <reason>". Renderer calls are identical; those three edits are the
# whole HTML difference (diffed against the previous code). The rule is checked
# by test_error_disclosure.py.
# 2026-10-06: every case with data moved, because method sentences now cite
# their method from `method_references`: the imputation line (here the
# numerical downshift, Tyanova et al. 2016) and the phosphosite Class I line
# (Olsen et al. 2006). Those citations are the whole HTML difference (diffed
# against the previous code for every mode); renderer calls are identical.
# 2026-10-08: every case that renders a table moved, because the report's
# tables now show display labels instead of the frames' column names (Mean,
# Lower whisker, ... instead of mean, lower_whisker, ...), round floats to two
# decimals, and the per-sample coverage table is titled "Coverage Table (per
# sample)" rather than "Coverage Summary Table". The labelling is covered by
# tests/repo/test_column_labels_in_sync.py; the empty-kinds cases are unchanged.
# Same day: the phosphosite coverage sections are titled "(per sample)" and
# "(per condition)" instead of "(Number)" and "(Quality)" (TOC and headings).
# 2026-10-09: the data-source line reads "Unfiltered Data" instead of "Data"
# (the bare word said nothing beside "Filtered Data"); every case with a
# provenance block moved.
# 2026-10-21: every case with data moved, because the volcano payload the report
# forwards gained `proteinCorrection` (false here; the phosphosite protein
# correction). Under this harness the probe's arguments are the
# figure's bytes, so the recorded `volcano_html` call and the HTML embedding it
# both move; stripping that one field from the payload reproduces every
# previous hash (verified for all 18 cases). The correction itself is covered
# by tests/unit/notebook/test_protein_correction.py and its report sentence by
# test_volcano_report_block.py.
# 2026-10-21, later: every case whose volcano renders interactively moved,
# because the report now forwards the entry's `highlight_outline` (true by
# default) to `volcano_html`, as it forwards the label counts to `volcano_png`;
# stripping that one keyword from the recorded call reproduces every previous
# hash (verified for all 18 cases). The option itself is covered by
# tests/unit/stats/test_volcano_highlight_options.py.
# 2026-10-30: Beta 2 changes default seeds to 1, retention-time defaults to
# scatter/100 bins, and phosphosite coverage to separate condition-count and
# localization plots. All 24 old fingerprints were reproduced against HEAD;
# HTML and renderer-call diffs contain only these intended changes.
SNAPSHOTS = {
    'qc_custom:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    'qc_custom:protein': 'a64884a1533f6944e066fb8bd3b297f1faffec8fb4d1a733fddc77384a068e38',
    'qc_custom:phospho': '367f4408e20d3dc1388436d7947af4a92edd4e6407ca44dbf4f7513f3c338fd2',
    'qc_custom:peptide,protein,phospho': '15d4acbeeec2e1e2fa4ee64012ff7c080219f4aac3ffe3e1b54c6c73502a3632',
    'qc_invalid:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    'qc_invalid:protein': '4d5cbd7c7d81985570ce3200fbd9217221b4d4ed6e92b56d3cc6ee0f655ec581',
    'qc_invalid:phospho': 'bcc7b12f6976395ce83ae502be54d8385097dac459b5494d3c8e335aaf345d73',
    'qc_invalid:peptide,protein,phospho': 'f2dae8ea2068f6b6fc1d722fb6a944d87bc7f6503566b8ee53caf3445c7bd55a',
    'default:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    # The three `default` snapshots below changed when advanced mode stopped
    # gating the appearance settings: the case's custom PCA title now reaches
    # the PCA renderer with the panel closed, as it reaches the plot on screen.
    # That forwarded title is the whole difference.
    'default:protein': '3a60785f8363b0b977f4d592b0ea2f5509f93175b3ceec9544fa097e440a8d65',
    'default:phospho': '951c73b044c230ffc3c43414ecac1b2d50a709c80e7ad77e7c99bf50d70f4e3a',
    'default:peptide,protein,phospho': 'c895ee6c8eb1a19613b4de6b4af58ca23a0bc64b593912b73b1cb536b1645ee7',
    'custom:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    'custom:protein': 'cc751016fcc8ca0759bf43a54e5e3683492be44afb76ebcec4d10b3a5ff0a2dc',
    'custom:phospho': 'dda30ee8923f2edcf1b45c77a2654b3196d7e9f61d200e646fcee9ad4ed3ca7b',
    'custom:peptide,protein,phospho': 'ffec78fd7625d4860089f575c988acd044aff351764bf3bb277ea27024da913c',
    'interactive:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    'interactive:protein': '9586f51b7f0ac8d4d542bbd24cebec03a47f3a93f7bfe76cb8e246a26f700346',
    'interactive:phospho': 'c6bea371ef65439ed82accf851a72da72c4cb9f921308c5ff8b924ba3cc7bddd',
    'interactive:peptide,protein,phospho': 'b01dd8a238fd451328f2405054613797373e0f316cd7beca3aa5c46ee0932b44',
    'failure:': '31e9b3e4680befa02db2ec0566567663337d540546ff762cdaa8127b85e1d920',
    'failure:protein': 'a7c5363a51694bdfb6dc138051f08da503cbbf8f347f850df1312febdd9abc5a',
    'failure:phospho': '6391509fe752427c4ba7e45e082e15aaf0270063241233c865e8e7ac5644ad77',
    'failure:peptide,protein,phospho': '6e7eebc844b194479e8c1af053d17ae479dbda1812f2d6363bb5cb94b58e3300',
}


def render_case(monkeypatch, kinds, mode, *, pca_scaling=None, selected_settings=False):
    monkeypatch.setattr(report, 'datetime', SimpleNamespace(now=lambda: datetime(2026, 9, 9)))
    comments = {key: [f'{key} <general> & text'] for key in report.GENERAL_COMMENTS}
    monkeypatch.setattr(report, 'GENERAL_COMMENTS', comments)
    monkeypatch.setattr(report, '_load_general_comments', lambda: comments)
    # 2026-09-30: the boxplot sections write their first sentences from the mode,
    # normalization and upload scale actually used; that wording is checked by
    # test_report_boxplot_text. The general comment is kept here so the hashes
    # keep proving the rest of the report is unchanged.
    monkeypatch.setattr(report, '_boxplot_comment_paragraphs',
                        lambda kind, mode, key: report.GENERAL_COMMENTS.get(key, []))
    annotation = SimpleNamespace(
        metadata=pd.DataFrame({'sample': ['A', 'B'], 'condition': ['A', 'B']}),
        source_data=pd.DataFrame({'value': [1, 2, 3]}), filtered_data=pd.DataFrame({'value': [1, 2]}),
        filter_config=SimpleNamespace(method='percentage', mode='overall', minPercentage=75, normalization='median'),
        imputation=SimpleNamespace(mode='numerical', downShift=1.8, width=0.3, seed=69, sampleWise=True),
    )
    monkeypatch.setattr(report, 'get_annotation', lambda kind: annotation if kind in kinds else None)
    monkeypatch.setattr(report, 'get_current_dataset', lambda kind: object() if kind in kinds else None)
    calls = []

    def renderer(name):
        def render(*args, **kwargs):
            values = [arg.model_dump() if hasattr(arg, 'model_dump') else arg for arg in args]
            encoded = json.dumps([name, values, kwargs], sort_keys=True)
            calls.append(encoded)
            if mode == 'failure' and 'pca' in name:
                raise ValueError('synthetic <failure> & details')
            if name.endswith('_table'):
                return pd.DataFrame({'Label': ['<A>&'], 'Value': [1.23456789]})
            return encoded if name.endswith('_html') or name == 'volcano_html' else encoded.encode()
        return render

    names = [name for name in vars(report) if name.startswith(('qc_', 'peptide_')) or
             name in ('completeness_missing_value_plot', 'completeness_missing_value_plot_peptide',
                      'phospho_coverage_png', 'phosphosite_plot_png', 'volcano_html')]
    for name in names:
        if callable(getattr(report, name)):
            monkeypatch.setattr(report, name, renderer(name))
    # The per-volcano counts come from run_volcano, which reads the real dataset
    # stores. Stub it deterministically so these fingerprints stay about report
    # structure; the counting itself is covered by test_volcano_report_block.py.
    monkeypatch.setattr(report, 'run_volcano', lambda payload: SimpleNamespace(
        upregulatedCount=312, downregulatedCount=288,
        totalRows=7450, notSignificantCount=6850))
    settings = {} if mode == 'default' else {
        'type': 'Interactive' if mode == 'interactive' else 'Normal', 'dataSource': 'imputed',
        'header': False, 'legend': False, 'plotDim': '3D', 'method': 'UMAP', 'addEllipses': True,
        'dotSize': 90, 'clusterMethod': 'KMeans', 'clusterCount': 1, 'colorBy': 'Sample',
        'umapNNeighbors': 999, 'umapMinDist': -1, 'tsnePerplexity': 2000, 'tsneLearningRate': 1,
        'dbscanEpsFactor': 0, 'dbscanMinSamples': 999, 'randomState': 123, 'widthCm': 90, 'heightCm': 2, 'dpi': 50,
    }
    if pca_scaling is not None:
        settings['scaling'] = pca_scaling
    context = SummaryReportContext(
        qcSettings={kind: {'pca': settings} for kind in ('protein', 'phospho')},
        advancedPlotMode=mode != 'default', qcAdvancedSettingsByTab={'pca': {'title': '<custom> & title'}},
        volcanoEntries=[SummaryVolcanoEntry(kind=kind, condition1='A', condition2='B', dataSource='imputed')
                        for kind in ('protein', 'phospho') if kind in kinds],
    )
    if mode in ('qc_custom', 'qc_invalid'):
        plot_settings = {
            'header': False, 'legend': False, 'widthCm': 32, 'heightCm': 18, 'dpi': 144,
            'mode': 'Median', 'outliers': True, 'includeId': True, 'text': True, 'textSize': 15,
            'method': 'Scatter', 'fullRange': True,
        }
        if mode == 'qc_invalid':
            plot_settings.update(header='false', legend='true', widthCm=-1, heightCm=100, dpi='invalid')
        for kind in ('protein', 'phospho'):
            context.qcSettings[kind].update({
                tab: plot_settings for tab in ('hist', 'box', 'cv', 'abundance', 'correlation')
            })
    if selected_settings:
        for kind in ('protein', 'phospho'):
            context.qcSettings[kind]['abundance'] = {'condition': 'B', 'type': 'Interactive' if mode == 'interactive' else 'Normal'}
        context.qcSettings['peptide'] = {'coverage': {'quantityColumn': 'LFQ intensity'}}
        context.completenessSettings = {'quantityColumn': 'LFQ intensity', 'binCount': 1, 'header': False,
            'text': False, 'textSize': 12, 'color': '#123456', 'widthCm': 25, 'heightCm': 15, 'dpi': 144}
    result = report.report_function(SummaryReportRequest(
        title='Report <title> & test', author='Author <name>', reportContext=context,
        textEntries={'Introduction': 'Intro <safe>\nsecond line', 'PCAProtAbove': 'Before <PCA>',
                     'PCAPhosBelow': 'After & PCA'},
    ))
    return result, calls


@pytest.mark.parametrize('kinds', [(), ('protein',), ('phospho',), ('peptide', 'protein', 'phospho')])
@pytest.mark.parametrize('mode', ['default', 'custom', 'interactive', 'failure', 'qc_custom', 'qc_invalid'])
def test_report_output_and_renderer_arguments_are_unchanged(monkeypatch, kinds, mode):
    monkeypatch.setattr(report, 'VERSION', 'Proteomics CoPYlot V1.0.0')
    result, calls = render_case(monkeypatch, kinds, mode)
    assert result.html.count(GENERATED_HTML_META) == 1
    # Only the newly added policy differs from the pre-hardening report baseline.
    result.html = result.html.replace('\n' + GENERATED_HTML_META, '', 1)
    # R7 intentionally changes only these sentences in this probe harness.
    # Their exact presence/meaning is checked by test_report_processing_provenance.
    # Removing them here proves all remaining HTML and renderer calls are unchanged.
    result.html = re.sub(
        r'<p>QC plots use the (?:normalized, )?unfiltered, non-imputed data unless another data source is stated. '
        r'Filtering and imputation apply to clustering and statistical plots according to '
        r'the selected data source; PCA also uses its stated data source.</p>', '', result.html
    )
    # 2026-09-24: normalization moved onto its own labelled line and a plot on
    # the unfiltered source now states the row counts instead of the filter
    # recipe. Both are wording checked by test_report_normalization_line; they
    # are removed here so the hashes keep proving the rest is unchanged.
    result.html = re.sub(r'<p><strong>Normalization:</strong>[^<]*</p>', '', result.html)
    result.html = re.sub(
        r'<p><strong>Filtering:</strong>[^<]*</p>', '<p><strong>Filtering:</strong></p>', result.html
    )
    # 2026-09-24: the PCA section states its parameters on a Projection line;
    # its wording is checked by test_report_pca_is_always_pca, so it is removed
    # here as the Normalization line is.
    result.html = re.sub(r'<p><strong>Projection:</strong>[^<]*</p>', '', result.html)
    # 2026-09-30: volcano blocks on imputed data state that imputed values were
    # tested as observations; the wording is checked by
    # test_imputed_value_counts, so it is removed here like the lines above.
    result.html = result.html.replace(
        '<p>' + html.escape(report.IMPUTED_AS_OBSERVED_SENTENCE) + '</p>', ''
    )
    # 2026-09-30: phosphosite sections state the session's localization
    # probability cutoff; the wording is checked by test_phospho_localization_cutoff.
    result.html = re.sub(r'<p><strong>Localization:</strong>[^<]*</p>', '', result.html)
    key = mode + ':' + ','.join(kinds)
    # Notebook availability is new response metadata, not a change to report contents.
    assert result.notebookToken is None and result.notebookError is None and result.notebookOmitted == []
    actual = json.dumps([result.model_dump(exclude={'notebookToken', 'notebookError', 'notebookOmitted'}), calls],
                        sort_keys=True).encode()
    assert hashlib.sha256(actual).hexdigest() == SNAPSHOTS[key]


@pytest.mark.parametrize('mode', ['default', 'interactive'])
@pytest.mark.parametrize('scaling', ['centered', 'unit_variance', 'auto'])
def test_report_pca_forwards_selected_scaling(monkeypatch, mode, scaling):
    result, calls = render_case(monkeypatch, ('protein', 'phospho'), mode, pca_scaling=scaling)
    projections = [json.loads(call) for call in calls if json.loads(call)[0] in ('qc_pca_plot', 'qc_pca_interactive_html')]
    assert len(projections) == 2
    assert all(call[2]['scaling'] == scaling for call in projections)
    assert all(call[0] == ('qc_pca_interactive_html' if mode == 'interactive' else 'qc_pca_plot') for call in projections)
    assert not result.warnings


@pytest.mark.parametrize('mode', ['default', 'interactive'])
def test_report_retains_abundance_quantity_and_missingness(monkeypatch, mode):
    result, calls = render_case(monkeypatch, ('peptide', 'protein', 'phospho'), mode, selected_settings=True)
    parsed = [json.loads(call) for call in calls]
    abundance = [call for call in parsed if call[0].startswith('qc_abundance')]
    assert len(abundance) == 2
    for name, _, options in abundance:
        assert options['condition'] == 'B'
        assert name == ('qc_abundance_interactive_html' if mode == 'interactive' else 'qc_abundance_plot')
        assert ('dpi' in options) == (mode != 'interactive')
    missing = [call for call in parsed if call[0].startswith('completeness_missing_value_plot')]
    assert len(missing) == 3
    for name, _, options in missing:
        assert options['bin_count'] == 1
        assert options['header'] is False and options['text'] is False
        assert options['text_size'] == 12 and options['color'] == '#123456'
        assert (options['width_cm'], options['height_cm'], options['dpi']) == (25, 15, 144)
        assert options.get('quantity_column') == ('LFQ intensity' if name.endswith('_peptide') else None)
    coverage = next(call for call in parsed if call[0] == 'qc_peptide_coverage_plot')
    assert coverage[2]['quantity_column'] == 'LFQ intensity'
    assert not result.warnings


def test_selected_quantity_changes_real_missingness_bars_and_stale_selection_warns(monkeypatch):
    from app.services import completeness_tools, peptide_tools
    frame = pd.DataFrame({'File.Name': ['r1', 'r2', 'r1', 'r2'],
        'Stripped.Sequence': ['AAA', 'AAA', 'BBB', 'BBB'],
        'Precursor.Quantity': [1, 1, 1, 1], 'LFQ intensity': [1, 0, 0, 1]})
    monkeypatch.setattr(peptide_tools, 'get_peptide_frame', lambda: frame)
    monkeypatch.setattr(completeness_tools, '_peptide_sample_columns', lambda _: ['r1', 'r2'])
    seen = []
    def capture(fig, plt, **kwargs):
        ax = fig.axes[0]
        seen.append(([label.get_text() for label in ax.get_xticklabels()], [bar.get_height() for bar in ax.patches]))
        plt.close(fig)
        return b'png'
    monkeypatch.setattr(completeness_tools, '_to_png_bytes', capture)
    report.completeness_missing_value_plot_peptide(**report._missing_plot_options({'quantityColumn': 'LFQ intensity'}, peptide=True))
    assert seen[-1] == (['1'], [2])
    report.completeness_missing_value_plot_peptide()
    assert seen[-1] == (['0'], [2])
    warnings = []
    report._safe_png(renderer=lambda: report.completeness_missing_value_plot_peptide(quantity_column='Removed intensity'),
                     fallback_title='Missingness', warnings=warnings)
    assert len(warnings) == 1 and 'Removed intensity' in warnings[0]
    assert len(seen) == 2


@pytest.mark.parametrize('render_type', ['Normal', 'Interactive'])
def test_abundance_stale_condition_warns_without_expanding_scope(monkeypatch, render_type):
    from app.services import plot_images
    monkeypatch.setattr(plot_images, '_abundance_mean_intensities', lambda _: (
        pd.DataFrame({'GeneNames': ['A'], 'B': [2.]}), ['B'], 'Protein', 'GeneNames'))
    warnings = []
    result = report._abundance_report('protein', {'condition': 'Removed', 'type': render_type}, None, warnings, 'Abundance')
    assert len(warnings) == 1 and 'Removed' in warnings[0]
    assert 'Could not render' in result


@pytest.mark.parametrize('settings', [{'binCount': 1000000000}, {'widthCm': float('inf')}, {'dpi': 100000}, {'binCount': 1.5}])
def test_missing_report_settings_reject_unsafe_render_sizes(settings):
    with pytest.raises(ValueError):
        report._missing_plot_options(settings)
