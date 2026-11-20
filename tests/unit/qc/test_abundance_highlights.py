import json

import numpy as np
import pandas as pd
import pytest

from app.services.abundance_highlights import (
    abundance_selection, abundance_interactive_highlights, abundance_static_highlights,
)


@pytest.fixture
def ranked():
    return pd.DataFrame({
        "Feature": ["AKT1_HUMAN;MTOR_HUMAN", "P~S1", "P~S2", "LOW", "MISSING"] * 2,
        "Condition": ["A"] * 5 + ["B"] * 5,
        "Rank": [1, 2, 3, 4, 5, 4, 3, 2, 1, 5],
        "log10Intensity": [4., 3., 2., 1., np.nan, 1., 2., 3., 4., np.inf],
    })


def test_selection_is_per_condition_and_does_not_change_values(ranked):
    original = ranked.copy(deep=True)
    assert ranked.index[abundance_selection(ranked, "mtor P~S1", 1, 1)].tolist() == [0, 1, 3, 5, 6, 8]
    assert ranked.index[abundance_selection(ranked, "P~S1")].tolist() == [1, 6]
    assert not abundance_selection(ranked, "P").any()
    assert abundance_selection(ranked, top_high=100, top_low=100).sum() == 8
    pd.testing.assert_frame_equal(ranked, original)


@pytest.mark.parametrize("outline,labels", [(True, False), (False, True), (True, True), (False, False)])
def test_interactive_overlays_preserve_condition_identity(ranked, outline, labels):
    import plotly.graph_objects as go
    fig = go.Figure()
    abundance_interactive_highlights(fig, ranked, {"A": "red", "B": "blue"}, {"A": "Control", "B": "Treatment"},
                                     "P~S1", 0, 0, outline, labels,
                                     {"abundanceHighlightLegendLabel": "Targets", "abundanceHighlightLegend": False})
    assert len(fig.data) == (2 if outline else 0)
    for trace, name, color in zip(fig.data, ["Control", "Treatment"], ["red", "blue"]):
        assert trace.legendgroup == name and trace.showlegend is False
        assert "text" not in trace.mode
        assert ("markers" in trace.mode) == outline
        assert trace.textfont.color == color
        assert trace.text is None


def test_static_labels_keep_both_conditions_and_one_label_per_point(ranked):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set(xlim=(0, 6), ylim=(0, 6))
    try:
        abundance_static_highlights(fig, ax, ranked, {"A": "red", "B": "blue"}, {"A": "Control", "B": "Treatment"},
                                    "LOW", 1, 1, True, True)
        texts = [text.get_text() for text in ax.texts]
        assert "LOW (Control)" in texts and "LOW (Treatment)" in texts
        assert len(texts) == len(set(texts)) == 4
        assert {text.get_color() for text in ax.texts} == {"red", "blue"}
        assert len(ax.collections[0].get_offsets()) == 4
    finally:
        plt.close(fig)


def test_outlines_follow_dot_sizes_and_use_a_separate_highlight_color(ranked):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    import plotly.graph_objects as go
    colors = {"A": "red", "B": "blue"}
    names = {"A": "Control", "B": "Treatment"}
    fig, ax = plt.subplots()
    try:
        ax.scatter([1], [1], s=4, label="Control")
        ax.scatter([1], [1], s=16, label="Treatment")
        abundance_static_highlights(fig, ax, ranked, colors, names, "P~S1", 0, 0, True, False)
        np.testing.assert_allclose(ax.collections[-1].get_sizes(), [4, 16])
        np.testing.assert_allclose(ax.collections[-1].get_edgecolors(), [to_rgba("#333333")] * 2)
        abundance_static_highlights(fig, ax, ranked, colors, names, "P~S1", 0, 0, True, False,
                                    {"abundanceHighlightWidth": 3, "abundanceHighlightColor": "black"})
        np.testing.assert_allclose(ax.collections[-1].get_sizes(), [4, 16])
        np.testing.assert_allclose(ax.collections[-1].get_edgecolors(), [to_rgba("black")] * 2)
    finally:
        plt.close(fig)
    figure = go.Figure()
    for name, size in [("Control", 4), ("Treatment", 9)]:
        figure.add_scatter(x=[1], y=[1], legendgroup=name, marker={"size": size})
    abundance_interactive_highlights(figure, ranked, colors, names, "P~S1", 0, 0, True, False)
    assert [trace.marker.size for trace in figure.data if trace.marker.symbol == "circle-open" and trace.legendgroup != "_copylot_highlight_legend"] == [4, 9]
    assert [trace.marker.color for trace in figure.data if trace.marker.symbol == "circle-open" and trace.legendgroup != "_copylot_highlight_legend"] == ["#333333", "#333333"]


@pytest.mark.parametrize("workflow,column", [("Protein", "ProteinNames"), ("Phosphosite", "PTM_Collapse_key"), ("Phosphoprotein", "Phosphoprotein")])
def test_notebook_options_preserve_ranks_and_use_gui_defaults(workflow, column):
    from proteomicscopylot import notebook_functions as nf
    frame = pd.DataFrame({column: ["HIGH", "MID", "LOW"], "a": [100., 10., 1.], "b": [1., 10., 100.]})
    if workflow == "Phosphoprotein":
        frame["PTM_Collapse_key"] = ["S1", "S2", "S3"]
    meta = pd.DataFrame({"sample": ["a", "b"], "condition": ["A", "B"]})
    default = nf.interactive_abundance_plot(frame, meta, "A", workflow, ["MID"], data_scale="linear")
    selected = [trace for trace in default.data if "Highlighted" in trace.name and trace.legendgroup != "_copylot_highlight_legend"]
    assert len(selected) == 1 and selected[0].mode == "markers"
    assert list(selected[0].x) == [2]
    assert selected[0].marker.symbol == "circle-open"
    assert selected[0].marker.color == "#333333"
    assert selected[0].marker.color != default.data[0].marker.color
    assert selected[0].marker.size == 6
    assert all("text" not in trace.mode for trace in default.data)
    explicit = nf.interactive_abundance_plot(frame, meta, "A", workflow, ["MID"], data_scale="linear",
                                             highlight_outline=True, highlight_labels=False)
    assert default.to_json() == explicit.to_json()
    figure = nf.interactive_abundance_plot(frame, meta, "All Conditions", workflow, data_scale="linear",
                                           top_high=1, top_low=1, search=["MID"], highlight_labels=True)
    overlays = [trace for trace in figure.data if "Highlighted" in trace.name and trace.legendgroup != "_copylot_highlight_legend"]
    assert len(overlays) == 2
    assert all(trace.text is None and trace.mode == "markers" for trace in overlays)
    assert all(list(trace.x) == [1, 2, 3] for trace in overlays)
    assert all(list(trace.y) == [2., 1., 0.] for trace in overlays)
    assert json.loads(figure.to_json())["data"]


def test_notebook_static_long_labels_stay_inside_axes():
    import matplotlib.pyplot as plt
    from proteomicscopylot import notebook_functions as nf
    frame = pd.DataFrame({"ProteinNames": ["H2A1B_HUMAN;H2A1D_HUMAN;H2A3_HUMAN", "LOW"],
                          "a": [100., 1.], "b": [90., 2.]})
    meta = pd.DataFrame({"sample": ["a", "b"], "condition": ["WT", "GRK"]})
    fig = nf.abundance_plot(frame, meta, data_scale="linear", top_high=1, top_low=1, highlight_labels=True)
    try:
        fig.canvas.draw()
        ax = fig.axes[0]
        assert len(ax.texts) >= 2
        for text in ax.texts:
            box = text.get_bbox_patch().get_window_extent(fig.canvas.get_renderer())
            assert ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
    finally:
        plt.close(fig)


def test_report_preserves_highlighting(monkeypatch):
    from app.services import report_function as report
    monkeypatch.setattr(report, "_safe_png", lambda renderer, **kwargs: renderer())
    monkeypatch.setattr(report, "report_plot", lambda renderer, **kwargs: kwargs)
    result = report._abundance_report("protein", {
        "search": "AKT1", "topHigh": 2, "topLow": 3,
        "highlightOutline": False, "highlightLabels": True,
    }, {"abundanceHighlightWidth": 2}, [], "Abundance")
    assert result["search"] == "AKT1"
    assert (result["top_high"], result["top_low"]) == (2, 3)
    assert result["highlight_outline"] is False and result["highlight_labels"] is True
    assert result["advanced"] == {"abundanceHighlightWidth": 2}

@pytest.mark.parametrize("term", ["GFRP", "GFRP_MOUSE", "gfrp"])
def test_mouse_suffix_matches_without_broad_phosphosite_matching(term):
    from proteomicscopylot import notebook_functions as nf
    frame = pd.DataFrame({"Feature": ["GFRP_MOUSE", "GFRP2_MOUSE", "GFRP_MOUSE~S1"],
                          "Condition": ["A"] * 3, "Rank": [1, 2, 3], "log10Intensity": [3., 2., 1.]})
    assert abundance_selection(frame, term).tolist() == [True, False, False]
    assert nf._abundance_selection(frame, term).tolist() == [True, False, False]


def test_gene_labels_preserve_gui_and_portable_ranks(monkeypatch):
    from app.services import plot_images
    from app.services.notebook_interactive import _abundance_module
    from app.services.notebook_qc import _qc_module
    data = pd.DataFrame({"ProteinNames": ["GFRP_MOUSE", "OTHER_MOUSE", "LAST_MOUSE"],
                         "GeneNames": ["Gfrp", "Gfrp", None], "a": [100., 10., 1.], "b": [1., 10., 100.]})
    meta = pd.DataFrame({"sample": ["a", "b"], "condition": ["A", "B"]})
    monkeypatch.setattr(plot_images, "linear_qc_frame", lambda kind: data)
    monkeypatch.setattr(plot_images, "_get_sample_columns", lambda kind, frame: ["a", "b"])
    monkeypatch.setattr(plot_images, "_metadata_for_kind", lambda *args: meta)
    original, _, _ = plot_images._abundance_rank_frame("protein")
    genes, _, _ = plot_images._abundance_rank_frame("protein", "genes")
    pd.testing.assert_frame_equal(genes[["Rank", "log10Intensity"]], original[["Rank", "log10Intensity"]])
    assert genes.Feature.tolist() == ["Gfrp", "Gfrp", "LAST_MOUSE"] * 2
    portable, _, _ = _abundance_module()._abundance_rank_frame(data, meta, "protein", "genes")
    pd.testing.assert_frame_equal(genes, portable)
    means, _, _, key = _qc_module("abundance")._abundance_mean_intensities(data, meta, "protein", "genes")
    assert key == "GeneNames" and len(means) == 3
    with pytest.raises(ValueError, match="Gene names"):
        _abundance_module()._abundance_rank_frame(data.drop(columns="GeneNames"), meta, "protein", "genes")


def test_abundance_labels_avoid_a_legend_inside_the_axes(ranked):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 6))
    try:
        ax.set(xlim=(0, 6), ylim=(0, 6))
        ax.scatter([1], [1], label="Control")
        legend = ax.legend(loc="center")
        abundance_static_highlights(fig, ax, ranked, {"A": "red", "B": "blue"}, {"A": "Control", "B": "Treatment"},
                                    "LOW", 1, 1, False, True)
        fig.canvas.draw()
        box = legend.get_window_extent(fig.canvas.get_renderer())
        assert ax.texts
        assert all(not text.get_bbox_patch().get_window_extent(fig.canvas.get_renderer()).overlaps(box) for text in ax.texts)
    finally:
        plt.close(fig)


def test_interactive_abundance_has_one_legend_key_and_retains_condition_groups(ranked):
    import plotly.graph_objects as go
    figure = go.Figure()
    for condition in ("A", "B"):
        figure.add_scatter(x=[1], y=[1], name=condition, legendgroup=condition, marker={"size": 6})
    abundance_interactive_highlights(figure, ranked, {"A": "red", "B": "blue"}, {"A": "A", "B": "B"},
                                     "LOW", 1, 1, True, False)
    entries = [trace for trace in figure.data if trace.showlegend and "Highlighted" in trace.name]
    assert [trace.name for trace in entries] == ["Highlighted"]
    assert list(entries[0].x) == [None]
    overlays = [trace for trace in figure.data if trace.name.endswith("(Highlighted)")]
    assert {trace.legendgroup for trace in overlays} == {"A", "B"}
    assert all(trace.showlegend is False for trace in overlays)
