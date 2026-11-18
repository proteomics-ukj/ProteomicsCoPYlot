import json
import subprocess
from types import SimpleNamespace

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd
import pytest

from app.services.condition_palette_store import get_category_palette, set_condition_palette
from app.services.copylot_core import COPYLOT_COLORS, MATPLOTLIB_TAB10, UKJ_COLORS, resolve_plot_colors
from app.services.notebook_source import bound_settings
from proteomicscopylot import notebook_functions as nb
from tests.repo_paths import ROOT


@pytest.fixture
def palette_conditions(monkeypatch):
    """Give each palette check its own active annotation, regardless of test order."""
    from app.services import annotation_store
    annotations = {}
    monkeypatch.setattr(annotation_store, 'get_annotation', lambda kind: annotations.get(kind))

    def set_conditions(kind, conditions):
        annotations[kind] = SimpleNamespace(metadata=pd.DataFrame({'condition': conditions}))

    return set_conditions


@pytest.mark.parametrize('name,colors', [('copylot', COPYLOT_COLORS), ('ukj', UKJ_COLORS), ('tab10', MATPLOTLIB_TAB10)])
def test_full_preset_is_retained_even_with_one_condition(name, colors, palette_conditions):
    from app.services.peptide_tools import peptide_rt_plot
    palette_conditions('protein', ['Only condition'])
    set_condition_palette('protein', {'Only condition': colors[0]})
    assert get_category_palette('protein') == list(colors)
    assert bound_settings(peptide_rt_plot)['plot_colors'] == list(colors)
    assert resolve_plot_colors(43, name) == [colors[i % len(colors)] for i in range(43)]
    assert bound_settings(peptide_rt_plot, plot_colors=['red'])['plot_colors'] == ['red']


def test_custom_category_palette_and_empty_palette(palette_conditions):
    palette_conditions('protein', ['A', 'B'])
    set_condition_palette('protein', {'A': '#123456', 'B': '#654321'})
    assert get_category_palette('protein') == ['#123456', '#654321']
    assert resolve_plot_colors(3, get_category_palette('protein')) == ['#123456', '#654321', '#123456']
    with pytest.raises(ValueError, match='At least one'):
        resolve_plot_colors(1, [])


@pytest.mark.parametrize('name,colors', [('copylot', COPYLOT_COLORS), ('ukj', UKJ_COLORS), ('tab10', MATPLOTLIB_TAB10)])
def test_notebook_points_keep_values_and_reference_lines(name, colors):
    data = pd.DataFrame({'A': [2., 3., 5., 8.], 'B': [3., 4., 7., 9.]})
    original = data.copy(deep=True)
    fig = nb.correlation_xy_intensity_single(data, 'A', 'B', plot_colors=name)
    try:
        assert to_hex(fig.axes[0].collections[0].get_facecolors()[0]).upper() == colors[0]
        assert fig.axes[0].lines[0].get_color() == 'red'
        np.testing.assert_array_equal(fig.axes[0].collections[0].get_offsets(), data[['A', 'B']])
        pd.testing.assert_frame_equal(data, original)
    finally:
        plt.close(fig)


def test_gui_sty_category_colors_do_not_shift_when_s_is_absent(monkeypatch, palette_conditions):
    from app.services import phospho_tools
    palette_conditions('phospho', ['Only condition'])
    set_condition_palette('phospho', {'Only condition': UKJ_COLORS[0]})
    monkeypatch.setattr(phospho_tools, 'phospho_sty_table', lambda: pd.DataFrame({
        'Residue': ['T', 'Y'], 'Count': [3, 2], 'Percentage': [60., 40.]}))
    seen = []

    def capture(fig, *args, **kwargs):
        seen.extend(to_hex(p.get_facecolor()).upper() for p in fig.axes[0].patches)
        plt.close(fig)
        return b'png'

    monkeypatch.setattr(phospho_tools, '_to_png_bytes', capture)
    phospho_tools.phospho_sty_png()
    assert seen == list(UKJ_COLORS[1:3])


def test_imputation_colors_do_not_change_seeded_results():
    data = pd.DataFrame({'A': [1., 3., np.nan, 7.], 'B': [2., np.nan, 4., 8.]})
    meta = pd.DataFrame({'sample': ['A', 'B'], 'condition': ['one', 'two']})
    expected = nb.impute_values(data, meta, seed=13)
    for palette in ('ukj', 'tab10'):
        pd.testing.assert_frame_equal(expected, nb.impute_values(data, meta, seed=13, plot_colors=palette))


@pytest.mark.parametrize('name', ['copylot', 'ukj', 'tab10'])
@pytest.mark.parametrize('plot', ['rt', 'qq'])
def test_gui_palette_keeps_reference_lines_and_marker_geometry(monkeypatch, name, plot, palette_conditions):
    from app.services import peptide_tools, plot_images
    palette_conditions('protein', ['A'])
    colors = resolve_plot_colors(3, name)
    set_condition_palette('protein', {'A': colors[0]})
    figures = []
    module = peptide_tools if plot == 'rt' else plot_images
    monkeypatch.setattr(module, '_to_png_bytes', lambda fig, *a, **kw: figures.append(fig) or b'png')
    if plot == 'rt':
        frame = pd.DataFrame({'Predicted.RT': [1., 2., 3.], 'RT': [1.1, 2.2, 2.9]})
        monkeypatch.setattr(module, 'get_peptide_frame', lambda: frame)
        module.peptide_rt_plot(add_line=True)
        expected_size, expected_alpha = 1., .2
    else:
        frame = pd.DataFrame({'A': [1., 2., 4., 8.]})
        monkeypatch.setattr(module, 'log2_qc_frame', lambda kind: frame)
        monkeypatch.setattr(module, '_get_sample_columns', lambda kind, data: ['A'])
        module.distribution_qqnorm_plot('protein')
        expected_size, expected_alpha = 2.5 ** 2, .5  # Matplotlib stores marker area.
    try:
        ax = figures[0].axes[0]
        dots = ax.collections[0]
        assert to_hex(dots.get_facecolors()[0]).upper() == colors[0]
        np.testing.assert_array_equal(dots.get_sizes(), [expected_size])
        assert dots.get_alpha() == expected_alpha
        assert len(ax.lines) == 1
        assert to_hex(ax.lines[0].get_color()) == '#ff0000'
        if plot == 'rt':
            np.testing.assert_array_equal(dots.get_offsets(), frame)
            np.testing.assert_array_equal(ax.lines[0].get_xdata(), [1., 3.])
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), [1., 3.])
    finally:
        for fig in figures:
            plt.close(fig)


@pytest.mark.requires_r('jsonlite')
def test_r_palettes_and_category_rendering_match_python(require_r, tmp_path):
    source = json.dumps((ROOT / 'proteomicscopylot/notebook_functions.R').as_posix())
    script = f'''source({source})
pdf(NULL)
for (name in c("copylot", "ukj", "tab10")) {{
  captured <- NULL
  barplot <- function(height, ...) {{
    args <- list(...); captured <<- args$col
    graphics::barplot(height, ...)
  }}
  data <- data.frame(PTM_Collapse_key=c("P1~G1_S1_M1", "P2~G2_T2_M1", "P3~G3_Y3_M1"),
    PTM_localization=c(.9,.8,.8))
  sty_plot(data, plot_colors=name)
  stopifnot(identical(unname(captured), .plot_colors(3,name)))
  cat(jsonlite::toJSON(list(name=name, colors=.plot_colors(43,name)), auto_unbox=TRUE), "\\n")
}}
dev.off()
'''
    path = tmp_path / 'palette.R'
    path.write_text(script, encoding='utf-8')
    result = subprocess.run([require_r('jsonlite'), '--vanilla', str(path)], capture_output=True, text=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    rows = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{')]
    assert len(rows) == 3
    for row in rows:
        assert [c.upper() for c in row['colors']] == resolve_plot_colors(43, row['name'])
