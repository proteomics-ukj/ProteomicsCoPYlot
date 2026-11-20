"""Top-N and typed terms share selection; rings follow dots and labels avoid legends."""
from __future__ import annotations

import json

import numpy as np
import pytest

from app.services import stats_tools


def _render(frame, **kwargs):
    defaults = dict(label_column="ProteinNames", condition1="A", condition2="B", p_value_threshold=0.05,
                    log2fc_threshold=1.0, use_uncorrected=False, highlight_terms=[])
    defaults.update(kwargs)
    return stats_tools._volcano_png_figure(frame, **defaults)


def _ring(ax):
    return [c for c in ax.collections if c.get_label() not in ("Upregulated", "Downregulated", "Not significant", "Unavailable")]


def _texts(fig):
    return {text.get_text() for text in fig.axes[0].texts if text.get_text()}


def test_outline_is_default_and_optional(volcano_frame) -> None:
    fig, plt = _render(volcano_frame, highlight_terms=["P0001", "P0025"])
    try:
        rings = _ring(fig.axes[0])
        assert len(rings) == 1 and len(rings[0].get_offsets()) == 2 and rings[0].get_label() == "Highlighted"
        assert _texts(fig) == set(), "a ring is not a label"
    finally:
        plt.close(fig)
    fig, plt = _render(volcano_frame, highlight_terms=["P0001", "P0025"], highlight_outline=False)
    try:
        assert _ring(fig.axes[0]) == [] and _texts(fig) == set()
    finally:
        plt.close(fig)


def test_labels_are_optional_and_union_with_the_top_hits(volcano_frame) -> None:
    top = set(stats_tools._volcano_top_hits(volcano_frame, label_column="ProteinNames", up=2, down=0)["ProteinNames"])
    highlighted = {"P0100", next(iter(top))}  # one non-hit and one that is already a top hit
    fig, plt = _render(volcano_frame, highlight_terms=sorted(highlighted), highlight_labels=True,
                       highlight_outline=False, label_top_up=2)
    try:
        drawn = _texts(fig)
        assert drawn <= top | highlighted and "P0100" in drawn
        assert len([t for t in fig.axes[0].texts if t.get_text() == next(iter(top))]) == 1, "labelled once"
        assert _ring(fig.axes[0]) == []
    finally:
        plt.close(fig)


def test_ring_size_follows_the_dot_and_width_does_not_change_diameter(volcano_frame) -> None:
    def ring_size(advanced):
        fig, plt = _render(volcano_frame, highlight_terms=["P0001"], advanced=advanced)
        try:
            return float(np.asarray(_ring(fig.axes[0])[0].get_sizes())[0])
        finally:
            plt.close(fig)
    default, general, own = ring_size(None), ring_size({"markerSize": 20}), ring_size({"markerSize": 20, "volcanoHighlightWidth": 3})
    assert general > default and own == general
    assert np.isclose(general, stats_tools.marker_area(20.0)) and np.isclose(own, stats_tools.marker_area(20.0))


def test_legend_entry_can_be_renamed_or_dropped(volcano_frame) -> None:
    fig, plt = _render(volcano_frame, highlight_terms=["P0001"],
                       advanced={"volcanoHighlightLegendLabel": "Kinases of interest"})
    try:
        assert "Kinases of interest" in [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    finally:
        plt.close(fig)
    fig, plt = _render(volcano_frame, highlight_terms=["P0001"], advanced={"volcanoHighlightLegend": False})
    try:
        assert "Highlighted" not in [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert len(_ring(fig.axes[0])) == 1, "the ring is still drawn"
    finally:
        plt.close(fig)


def test_interactive_plot_honours_outline_and_legend_but_never_labels(volcano_frame) -> None:
    frame = volcano_frame.assign(pval=0.01, adj_pval=0.02, selected_pval=0.02)
    figure = stats_tools._volcano_figure(frame, label_column="ProteinNames", condition1="A", condition2="B",
                                         p_value_threshold=0.05, log2fc_threshold=1.0, use_uncorrected=False,
                                         highlight_terms=["P0001"], highlight_labels=True,
                                         advanced={"volcanoHighlightLegendLabel": "Mine", "volcanoHighlightLegend": False})
    traces = json.loads(figure.to_json())["data"]
    highlighted = [t for t in traces if t["name"].endswith("(Mine)")]
    assert len(highlighted) == 1 and highlighted[0].get("showlegend") is False
    assert not any(t.get("mode") == "text" or "text" in str(t.get("mode", "")) for t in traces if t["name"].endswith("(Mine)"))
    figure = stats_tools._volcano_figure(frame, label_column="ProteinNames", condition1="A", condition2="B",
                                         p_value_threshold=0.05, log2fc_threshold=1.0, use_uncorrected=False,
                                         highlight_terms=["P0001"], highlight_outline=False)
    assert not any("(" in t["name"] for t in json.loads(figure.to_json())["data"])


@pytest.mark.parametrize("outline,labels", [(True, True), (False, False)])
def test_the_report_entry_carries_both_switches(outline, labels) -> None:
    from app.schemas.summary import SummaryVolcanoEntry
    entry = SummaryVolcanoEntry(kind="protein", condition1="A", condition2="B", highlightTerms=["X"],
                                highlightOutline=outline, highlightLabels=labels)
    assert (entry.highlightOutline, entry.highlightLabels) == (outline, labels)
    assert SummaryVolcanoEntry(kind="protein", condition1="A", condition2="B").highlightOutline is True


def test_every_volcano_route_takes_the_switches_the_panel_sends() -> None:
    """The panel sends highlightOutline to all four routes and highlightLabels to the
    static ones; a route without the parameter silently ignores it (found on the plain
    interactive route, 2026-10-22)."""
    from fastapi.testclient import TestClient
    from app.main import app
    paths = TestClient(app).get("/openapi.json").json()["paths"]
    for route, static in (("volcano.html", False), ("volcano-control.html", False),
                          ("volcano.png", True), ("volcano-control.png", True)):
        names = {p["name"] for p in paths[f"/api/plots/stats/{{kind}}/{route}"]["get"]["parameters"]}
        assert "highlightOutline" in names, route
        assert ("highlightLabels" in names) is static, route


def test_top_hits_and_custom_terms_share_outline_and_label_switches(volcano_frame):
    top = stats_tools._volcano_top_hits(volcano_frame, label_column="ProteinNames", up=2, down=2)
    expected = set(top.index) | set(volcano_frame.index[volcano_frame.ProteinNames == "P0100"])
    for outline in (True, False):
        for labels in (True, False):
            fig, plt = _render(volcano_frame, label_top_up=2, label_top_down=2, highlight_terms=["P0100"],
                              highlight_outline=outline, highlight_labels=labels)
            try:
                assert sum(len(ring.get_offsets()) for ring in _ring(fig.axes[0])) == (len(expected) if outline else 0)
                assert bool(_texts(fig)) == labels
                legend_box = fig.axes[0].get_legend().get_window_extent(fig.canvas.get_renderer())
                assert all(not text.get_window_extent(fig.canvas.get_renderer()).overlaps(legend_box) for text in fig.axes[0].texts)
            finally:
                plt.close(fig)
        interactive = stats_tools._volcano_figure(volcano_frame, label_column="ProteinNames", condition1="A", condition2="B",
            p_value_threshold=.05, log2fc_threshold=1, use_uncorrected=False, highlight_terms=["P0100"],
            label_top_up=2, label_top_down=2, highlight_outline=outline, advanced={"markerSize": 13, "volcanoHighlightWidth": 2.5})
        selected = [trace for trace in interactive.data if "Highlighted" in trace.name and trace.legendgroup != "_copylot_highlight_legend"]
        assert sum(len(trace.x) for trace in selected) == (len(expected) if outline else 0)
        assert all(trace.marker.size == 13 and trace.marker.line.width == 2.5 for trace in selected)


def test_mouse_highlight_normalization_preserves_sites():
    assert stats_tools._normalize_highlight_token("GFRP_MOUSE") == stats_tools._normalize_highlight_token("GFRP")
    assert stats_tools._normalize_highlight_token("GFRP_MOUSE~S1") != stats_tools._normalize_highlight_token("GFRP")


def test_volcano_labels_stay_within_both_axes_and_one_interactive_legend_key(volcano_frame):
    fig, plt = _render(volcano_frame, label_top_up=20, label_top_down=20, highlight_labels=True)
    try:
        fig.canvas.draw()
        ax = fig.axes[0]
        assert ax.texts
        for text in ax.texts:
            box = text.get_window_extent(fig.canvas.get_renderer())
            assert ax.bbox.contains(box.x0, box.y0) and ax.bbox.contains(box.x1, box.y1)
    finally:
        plt.close(fig)
    figure = stats_tools._volcano_figure(volcano_frame, label_column="ProteinNames", condition1="A", condition2="B",
        p_value_threshold=.05, log2fc_threshold=1, use_uncorrected=False, highlight_terms=["P0100"],
        label_top_up=2, label_top_down=2)
    assert [trace.name for trace in figure.data if trace.showlegend and "Highlighted" in trace.name] == ["Highlighted"]
    for trace in figure.data:
        if trace.name.endswith("(Highlighted)"):
            assert trace.showlegend is False
            assert trace.legendgroup in {"Upregulated", "Downregulated", "Not significant"}
            assert any(base.name == trace.legendgroup and base.showlegend for base in figure.data)
