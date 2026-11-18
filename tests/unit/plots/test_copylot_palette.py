"""The approved palette must agree across the app and standalone notebooks."""
import hashlib
import json
import re
import subprocess

import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
import pandas as pd
import pytest

from app.services.condition_palette_store import build_condition_color_map, set_condition_palette
from app.services.r_runtime import detect_r_runtime
from proteomicscopylot._core import COPYLOT_COLORS, copylot_colors
from proteomicscopylot import notebook_functions as notebook
from tests.repo_paths import ROOT


def test_approved_palette_is_fixed_and_repeats_after_forty():
    assert len(COPYLOT_COLORS) == len(set(COPYLOT_COLORS)) == 40
    assert copylot_colors(0) == []
    assert copylot_colors(81) == list(COPYLOT_COLORS) * 2 + [COPYLOT_COLORS[0]]
    # Pin every approved color, not only the brand opening.
    assert hashlib.sha256(''.join(COPYLOT_COLORS).encode()).hexdigest() == '616b02967bce21d14cb8e5186bdfa01d5c6d2d1937339b2e48d301687baa8c33'
    ts = (ROOT / 'frontend/src/lib/conditionColors.ts').read_text()
    assert re.findall(r'"(#[0-9A-F]{6})"', ts) == list(COPYLOT_COLORS)


def test_gui_palette_cycles_but_explicit_override_wins():
    names = [f'C{i}' for i in range(81)]
    set_condition_palette('protein', {'C40': '#123456'})
    try:
        colors = build_condition_color_map('protein', names)
        expected = copylot_colors(81)
        expected[40] = '#123456'
        assert list(colors.values()) == expected
    finally:
        set_condition_palette('protein', {})


@pytest.mark.requires_r()
def test_r_palette_matches_python_through_multiple_cycles():
    path = json.dumps((ROOT / 'proteomicscopylot/notebook_functions.R').as_posix())
    result = subprocess.run([detect_r_runtime()['executable'], '--vanilla', '-e',
        f'source({path}); cat(paste(.plot_colors(81),collapse=",")); stopifnot(identical(.plot_colors(3,c("red","blue")),c("red","blue","red")))'],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout.strip().split(',') == copylot_colors(81)


def test_abundance_static_and_interactive_keep_all_41_conditions_and_values():
    samples = [f'S{i}' for i in range(41)]
    meta = pd.DataFrame({'sample': samples, 'condition': samples})
    data = pd.DataFrame(np.random.default_rng(13).uniform(10,100,(6,41)),columns=samples)
    data.insert(0,'ProteinNames',[f'P{i}' for i in range(6)])
    original=data.copy(deep=True)
    figures=[]
    try:
        for custom in (None,['#123456']):
            fig=notebook.abundance_plot(data,meta,data_scale='linear',plot_colors=custom,highlight_outline=False,height_cm=36)
            figures.append(fig)
            static=[to_hex(c.get_facecolors()[0]).upper() for c in fig.axes[0].collections]
            interactive=notebook.interactive_abundance_plot(data,meta,'All Conditions',data_scale='linear',plot_colors=custom,highlight_outline=False)
            traces_by_condition={t.name: t for t in interactive.data if t.mode=='markers'}
            traces=[traces_by_condition[name] for name in samples]
            expected=copylot_colors(41) if custom is None else custom*41
            assert static == expected
            assert [t.marker.color for t in traces] == expected
            for collection,trace in zip(fig.axes[0].collections,traces):
                points = collection.get_offsets()
                points = points[np.argsort(points[:, 0])]
                np.testing.assert_allclose(points[:, 0], trace.x)
                np.testing.assert_allclose(points[:, 1], trace.y)
        pd.testing.assert_frame_equal(data,original)
    finally:
        for fig in figures: plt.close(fig)
