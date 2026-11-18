"""Behavioral counterparts for the intentionally different private helper APIs."""
import json
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from app.services.abundance_highlights import abundance_selection
from proteomicscopylot import notebook_functions as notebook
from proteomicscopylot._core import copylot_colors
from tests.repo_paths import ROOT


def r_json(require_r, script, *packages):
    source = json.dumps((ROOT / 'proteomicscopylot/notebook_functions.R').as_posix())
    result = subprocess.run(
        [require_r('jsonlite', *packages), '--vanilla', '-'],
        input=f'source({source})\n{script}', text=True, capture_output=True, timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return json.loads(result.stdout)


@pytest.mark.requires_r('jsonlite')
def test_copylot_color_function_matches_r_vector_repetition(require_r):
    counts = [0, 1, 2, 40, 41, 80]
    actual = r_json(require_r, 'cat(jsonlite::toJSON(lapply(c(0,1,2,40,41,80), .plot_colors)))')
    assert actual == [copylot_colors(n) for n in counts]


@pytest.mark.requires_r('jsonlite')
def test_abundance_selection_matches_across_three_pipelines(require_r):
    cases = [
        {'search': search, 'high': high, 'low': low}
        for search in (None, '', 'mtor P~S1', ['mtor', 'P~S1'], ['mtor P~S1'],
                       [' MTOR_human\tP~S1 ', 'LOW'], 'P', 'P~S1')
        for high, low in ((0, 0), (1, 1), (100, 100))
    ]
    frame = pd.DataFrame({
        'Feature': ['AKT1_HUMAN;MTOR_HUMAN', 'P~S1', 'P~S2', 'LOW', 'MISSING', 'INF'] * 2,
        'Condition': ['A'] * 6 + ['B'] * 6,
        'log10Intensity': [4., 3., 3., 1., np.nan, np.inf, 1., 3., 3., 4., np.nan, -np.inf],
    })
    frame['Rank'] = frame.groupby('Condition')['log10Intensity'].rank(ascending=False, method='first')
    original = frame.copy(deep=True)
    script = f'''
cases <- jsonlite::fromJSON({json.dumps(json.dumps(cases))}, simplifyVector=FALSE)
features <- c("AKT1_HUMAN;MTOR_HUMAN","P~S1","P~S2","LOW","MISSING","INF")
result <- lapply(cases, function(x) c(
  .abundance_selection(features,c(4,3,3,1,NA,Inf),unlist(x$search),x$high,x$low),
  .abundance_selection(features,c(1,3,3,4,NA,-Inf),unlist(x$search),x$high,x$low)))
cat(jsonlite::toJSON(result))
'''
    actual = r_json(require_r, script)
    for case, r_mask in zip(cases, actual, strict=True):
        options = dict(search=case['search'], top_high=case['high'], top_low=case['low'])
        gui = abundance_selection(frame, **options)
        python = notebook._abundance_selection(frame, **options)
        assert gui.tolist() == python.tolist() == r_mask, case
    # Assert the intended match as well as agreement between implementations.
    assert frame.index[notebook._abundance_selection(frame, ['mtor P~S1'])].tolist() == [0, 1, 6, 7]
    pd.testing.assert_frame_equal(frame, original)


@pytest.mark.requires_r('jsonlite', 'plotly')
def test_tied_abundances_keep_input_order_in_static_and_interactive_plots(require_r, monkeypatch):
    values = [4,3,2,1,3,1,3,4,4,3,1,3,1,2,1,4,1,4,2,2,2,1,4,4,4,2,4,3,3,3]
    data = pd.DataFrame({'ProteinNames': [f'P{i:02}' for i in range(30)], 'a': values, 'b': values[::-1]})
    meta = pd.DataFrame({'sample': ['a', 'b'], 'condition': ['A', 'B']})
    original = data.copy(deep=True)
    script = f'''
v <- c({','.join(map(str, values))})
d <- data.frame(ProteinNames=sprintf("P%02d",0:29),a=v,b=rev(v))
m <- data.frame(sample=c("a","b"),condition=c("A","B"))
f <- plotly::plotly_build(interactive_abundance_plot(d,m,"All Conditions",data_scale="linear",
  top_high=1,top_low=1,highlight_labels=TRUE))
traces <- Filter(function(t) grepl("Highlighted",t$name) && t$legendgroup != "_copylot_highlight_legend",f$x$data)
cat(jsonlite::toJSON(lapply(traces,function(t) list(labels=sub("<br>.*$", "", t$hovertext),x=t$x,y=t$y))))
'''
    r_traces = r_json(require_r, script, 'plotly')
    figure = notebook.interactive_abundance_plot(data, meta, 'All Conditions', data_scale='linear',
                                               top_high=1, top_low=1, highlight_labels=True)
    overlays = [trace for trace in figure.data if 'Highlighted' in trace.name and trace.legendgroup != '_copylot_highlight_legend']
    assert len(overlays) == len(r_traces) == 2
    expected_labels = []
    for sample, condition, trace, r_trace in zip(['a', 'b'], ['A', 'B'], overlays, r_traces, strict=True):
        ordered = data.sort_values(sample, ascending=False, kind='stable')
        selected = ordered.iloc[[0, -1]]
        labels = [f'{feature} ({condition})' for feature in selected.ProteinNames]
        expected_labels.extend(labels)
        assert trace.text is None
        assert [row[0] for row in trace.customdata] == r_trace['labels'] == selected.ProteinNames.tolist()
        assert list(trace.x) == r_trace['x'] == [1, 30]
        np.testing.assert_allclose(trace.y, np.log10(selected[sample]))
        np.testing.assert_allclose(trace.y, r_trace['y'], rtol=1e-4)
        # Base points retain the same abundance values, now ordered stably for ties.
        base = next(t for t in figure.data if t.name == condition)
        np.testing.assert_allclose(base.y, np.log10(ordered[sample]))
    assert expected_labels[:2] == ['P00 (A)', 'P21 (A)']
    # Exercise the GUI's single-condition static path, not just its selection helper.
    from app.services import plot_images
    means = pd.DataFrame({'ProteinNames': data.ProteinNames, 'A': np.log10(values)})
    monkeypatch.setattr(plot_images, '_abundance_mean_intensities', lambda *a, **k: (means, ['A'], 'Protein', 'ProteinNames'))
    monkeypatch.setattr(plot_images, '_to_png_bytes', lambda fig, *a, **k: fig)
    gui = plot_images.qc_abundance_plot('protein', condition='A', top_high=1, top_low=1, highlight_labels=True)
    try:
        assert sorted(t.get_text() for t in gui.axes[0].texts) == ['P00', 'P21']
    finally:
        plt.close(gui)
    static = notebook.abundance_plot(data, meta, data_scale='linear', top_high=1, top_low=1,
                                    highlight_labels=True)
    try:
        assert sorted(t.get_text() for t in static.axes[0].texts) == sorted(expected_labels)
    finally:
        plt.close(static)
    pd.testing.assert_frame_equal(data, original)


@pytest.mark.requires_r('jsonlite', 'plotly')
@pytest.mark.parametrize('features', [
    ['GENE_1', 'GENE_2', 'GENE_OTHER', 'GFRP_MOUSE', 'GFRP_MOUSE~S1'],
    ['P12345~S10_M1', 'P12345~S10_M2', 'NP_1234', 'NP_5678', 'OTHER'],
])
def test_distinct_identifier_suffixes_are_not_highlighted_together(require_r, features):
    from app.services import stats_tools
    frame = pd.DataFrame({'Feature': features, 'Condition': 'A', 'Rank': range(1, 6),
                          'log10Intensity': [5., 4., 3., 2., 1.]})
    script = f'features <- jsonlite::fromJSON({json.dumps(json.dumps(features))})\n' + '''
v <- data.frame(ProteinNames=features,log2FC=1:5,neg_log10_selected_pval=5:1,significance="Upregulated")
f <- plotly::plotly_build(volcano_plot_h(v,"A","B",highlight_terms=features[1]))
traces <- Filter(function(t) grepl("(highlighted)",t$name,fixed=TRUE),f$x$data)
cat(jsonlite::toJSON(list(mask=.abundance_selection(features,5:1,features[1]),x=unlist(lapply(traces,function(t)t$x)))))
'''
    actual = r_json(require_r, script, 'plotly')
    expected = [True, False, False, False, False]
    assert abundance_selection(frame, features[0]).tolist() == expected
    assert notebook._abundance_selection(frame, features[0]).tolist() == expected
    assert actual['mask'] == expected
    volcano = pd.DataFrame({'ProteinNames': features, 'log2FC': range(1, 6),
                            'neg_log10_selected_pval': [5., 4., 3., 2., 1.], 'significance': 'Upregulated'})
    figure = notebook.volcano_plot_h(volcano, 'A', 'B', highlight_terms=[features[0]])
    assert [x for t in figure.data if '(highlighted)' in t.name for x in t.x] == actual['x'] == [1]
    gui, _ = stats_tools._volcano_png_figure(volcano, label_column='ProteinNames', condition1='A', condition2='B',
        p_value_threshold=.05, log2fc_threshold=1, use_uncorrected=False, highlight_terms=[features[0]])
    try:
        rings = [c for c in gui.axes[0].collections if c.get_label() == 'Highlighted']
        assert rings[0].get_offsets()[:, 0].tolist() == [1]
    finally:
        plt.close(gui)


@pytest.mark.requires_r('jsonlite', 'plotly')
@pytest.mark.parametrize('palette', ['copylot', 'ukj', 'tab10', ['#123456']])
def test_interactive_abundance_resolves_named_and_repeated_palettes(require_r, palette, monkeypatch):
    from proteomicscopylot._core import resolve_plot_colors
    data = pd.DataFrame({'ProteinNames': ['P1', 'P2'], 'a': [10, 20], 'b': [30, 40]})
    meta = pd.DataFrame({'sample': ['a', 'b'], 'condition': ['A', 'B']})
    script = f'''
d <- data.frame(ProteinNames=c("P1","P2"),a=c(10,20),b=c(30,40))
m <- data.frame(sample=c("a","b"),condition=c("A","B"))
p <- jsonlite::fromJSON({json.dumps(json.dumps(palette))})
f <- plotly::plotly_build(interactive_abundance_plot(d,m,"All Conditions",data_scale="linear",plot_colors=p))
cat(jsonlite::toJSON(setNames(lapply(f$x$data,function(t)toupper(t$marker$color)),vapply(f$x$data,function(t)t$name,character(1))),auto_unbox=TRUE))
'''
    actual = r_json(require_r, script, 'plotly')
    figure = notebook.interactive_abundance_plot(data, meta, 'All Conditions', data_scale='linear', plot_colors=palette)
    expected = dict(zip(['A', 'B'], resolve_plot_colors(2, palette)))
    # Plotly may order traces by first occurrence after intensity sorting;
    # compare the condition/color association, not trace order or hex case.
    assert {t.name: t.marker.color.upper() for t in figure.data} == actual == expected
    from app.services import plot_images
    from app.services.condition_palette_store import set_condition_palette
    set_condition_palette('protein', expected)
    ranks = pd.DataFrame({'Feature': ['P1', 'P2'] * 2, 'Condition': ['A', 'A', 'B', 'B'],
                          'Rank': [2, 1, 2, 1], 'log10Intensity': np.log10([10, 20, 30, 40])})
    monkeypatch.setattr(plot_images, '_abundance_rank_frame', lambda *a, **k: (ranks, ['A', 'B'], 'Protein'))
    monkeypatch.setattr(plot_images, '_plotly_html', lambda fig: fig)
    gui = plot_images.qc_abundance_interactive_html('protein')
    assert {t.name: t.marker.color.upper() for t in gui.data} == expected
