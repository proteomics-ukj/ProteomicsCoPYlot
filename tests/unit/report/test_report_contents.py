"""The report's contents, and the plots a user added to it.

The contents are written from the headings, so the two cannot disagree. Added
plots form their own sections behind the fixed ones, are drawn the way they
were added (advanced appearance and size included), and a plot that cannot be
drawn is left out with one warning.
"""

from __future__ import annotations

import html as html_module
import re

import pytest

import app.services.report_function as report
from app.schemas.summary import SummaryPlotEntry, SummaryReportContext, SummaryReportRequest, SummaryVolcanoEntry
from app.services.annotation_store import get_annotation

pytestmark = pytest.mark.usefixtures("example_project")

VOLCANO = {"dataSource": "data", "condition1": "WT", "condition2": "GRK", "pValueThreshold": 0.05,
           "log2fcThreshold": 1.0, "testType": "unpaired", "useUncorrected": False}
RENDER = {"type": "Normal", "widthCm": 11, "heightCm": 7, "dpi": 96}


def _entry(entry_type: str, request: dict, *, kind: str = "protein", render: dict | None = None,
           advanced: dict | None = None) -> SummaryPlotEntry:
    return SummaryPlotEntry(type=entry_type, kind=kind, request=request, render=render or RENDER, advanced=advanced)


@pytest.fixture
def figures(monkeypatch):
    """Stand in for every figure and record how each was asked for."""
    calls: list[tuple[str, tuple, dict]] = []

    def draw(renderer, *args, **options):
        calls.append((renderer.__name__, args, options))
        return "<div>figure</div>" if renderer.__name__.endswith("_html") else b"png"

    monkeypatch.setattr(report, "report_plot", draw)
    return calls


def _build(entries=(), **context):
    return report.report_function(SummaryReportRequest(
        title="T", author="A", reportContext=SummaryReportContext(entries=list(entries), **context)))


def _contents(markup: str) -> list[tuple[str, str]]:
    return [(anchor, html_module.unescape(text))
            for anchor, text in re.findall(r"<p class='toc-(?:group|item)'><a href='#([^']+)'>([^<]*)</a></p>", markup)]


def _headings(markup: str) -> list[tuple[str, str]]:
    return [(anchor, html_module.unescape(text))
            for anchor, text in re.findall(r"<h[12][^>]* id='([^']+)'>([^<]*)</h[12]>", markup)]


def _samples() -> list[str]:
    return get_annotation("protein").metadata["sample"].astype(str).tolist()


def _all_entries() -> list[SummaryPlotEntry]:
    first, second = _samples()[:2]
    protein = str(get_annotation("protein").source_data["ProteinNames"].iloc[0])
    return [
        _entry("single_line", {"proteins": [protein], "conditions": ["WT", "GRK"], "dataSource": "data"}),
        _entry("single_box", {"protein": protein, "conditions": ["WT", "GRK"], "dataSource": "data"}),
        _entry("anova", {"source": "anova", "dataSource": "data", "pValueThreshold": 0.05}),
        _entry("volcano", VOLCANO),
        _entry("venn", {"mode": "Condition", "first": "WT", "second": "GRK", "minPresent": 2}),
        _entry("pearson", {"mode": "Single", "sample1": first, "sample2": second}),
        _entry("volcano", VOLCANO, kind="phospho"),
    ]


def test_contents_are_the_headings_and_every_line_is_a_link(figures) -> None:
    for entries in ([], _all_entries()):
        markup = _build(entries).html
        assert _contents(markup) == _headings(markup)
        assert len(_contents(markup)) > 10
        anchors = [anchor for anchor, _ in _headings(markup)]
        assert len(anchors) == len(set(anchors))


def test_added_plots_follow_the_fixed_sections_in_report_order(figures) -> None:
    titles = [text for _, text in _contents(_build(_all_entries()).html)]
    protein = titles[titles.index("1 Protein Level") + 1:titles.index("2 Phosphosite Level")]
    assert protein[8:] == ["1.9 Pearson Correlation", "1.10 Venn Diagram", "1.11 Volcano Plot",
                           "1.12 One-way ANOVA", "1.13 Single Protein"]
    phospho = titles[titles.index("2 Phosphosite Level") + 1:]
    assert phospho[-2:] == ["2.11 Correlation Plot", "2.12 Volcano Plot"]


def test_boxplots_come_before_line_plots_and_every_plot_states_its_data(figures) -> None:
    markup = _build(_all_entries()).html
    section = markup[markup.index("1.13 Single Protein</h2>"):markup.index("<h1 class='level' id='s-2'>")]
    assert section.index("<h3>Boxplot 1:") < section.index("<h3>Line plot 1:")
    assert section.count("<strong>Data source:</strong> Unfiltered Data") == 2
    for name in ("Pearson", "Venn", "Anova", "Single"):
        assert report.GENERAL_COMMENTS[f"{name}Prot"] and report.GENERAL_COMMENTS[f"{name}Phos"]
    assert report.GENERAL_COMMENTS["OraProt"]


def test_an_entry_is_drawn_with_its_own_appearance_and_size(figures) -> None:
    first, second = _samples()[:2]
    _build(
        [_entry("volcano", VOLCANO, advanced={"title": "Mine"}),
         _entry("pearson", {"mode": "Single", "sample1": first, "sample2": second},
                render={"type": "Interactive", "widthCm": 11, "heightCm": 7, "dpi": 96}, advanced={"title": "Scatter"})],
        advancedPlotSettingsByKey={"stats:volcano": {"title": "Live setting"}},
    )
    volcano = next(options for name, _, options in figures if name == "volcano_png")
    assert volcano["advanced"] == {"title": "Mine"}
    assert (volcano["width_cm"], volcano["height_cm"], volcano["dpi"]) == (11, 7, 96)
    pearson = next(options for name, _, options in figures if name == "pearson_correlation_html")
    assert pearson["advanced"] == {"title": "Scatter"} and "dpi" not in pearson


def test_volcano_entries_of_an_older_browser_tab_keep_the_live_appearance(figures) -> None:
    result = _build(
        volcanoEntries=[SummaryVolcanoEntry(kind="protein", condition1="WT", condition2="GRK", renderType="Normal")],
        advancedPlotSettingsByKey={"stats:volcano": {"title": "Live setting"}},
    )
    assert "1.9 Volcano Plot</h2>" in result.html
    volcano = next(options for name, _, options in figures if name == "volcano_png")
    assert volcano["advanced"] == {"title": "Live setting"} and "width_cm" not in volcano


def test_a_plot_that_cannot_be_drawn_is_left_out_with_one_warning(figures) -> None:
    result = _build([
        _entry("volcano", {**VOLCANO, "condition2": "WT"}),
        _entry("venn", {"mode": "Condition", "first": "WT", "second": "Removed"}),
        _entry("volcano", VOLCANO),
    ])
    assert len(result.warnings) == 2, result.warnings
    assert all("could not be rendered" in warning for warning in result.warnings)
    titles = [text for _, text in _contents(result.html)]
    assert "1.9 Volcano Plot" in titles and not any("Venn" in title for title in titles)
    assert result.html.count("<h3 id='e-protein-volcano-") == 1
    assert "Could not render" not in result.html


def test_an_ora_links_to_its_volcano_only_when_that_is_in_the_report(figures) -> None:
    ora = {**VOLCANO, "source": "volcano", "direction": "up", "adjPValueCutoff": 0.05, "topN": 5,
           "minTermSize": 20, "maxTermSize": 300, "backgroundMode": "experiment", "organism": "human"}
    alone = _build([_entry("ora", ora)]).html
    assert "GO Enrichment (ORA)</h2>" in alone and "Go to the underlying plot" not in alone
    assert "<strong>Source:</strong> upregulated in the comparison WT vs GRK" in alone
    linked = _build([_entry("ora", ora), _entry("volcano", {**VOLCANO, "highlightTerms": ["P1"]})]).html
    target = re.search(r"<a href='#(e-protein-volcano-\d+)'>Go to the underlying plot.</a>", linked)
    assert target and f"<h3 id='{target.group(1)}'>" in linked
    other = _build([_entry("ora", ora), _entry("volcano", {**VOLCANO, "log2fcThreshold": 2.0})]).html
    assert "Go to the underlying plot" not in other


def test_the_pca_follows_its_tab(figures, monkeypatch) -> None:
    monkeypatch.setattr(report, "most_processed_data_source", lambda kind: "filtered")
    _build()
    assert all(options["data_source"] == "filtered" for name, _, options in figures if name == "qc_pca_plot")
    figures.clear()
    _build(qcSettings={"protein": {"pca": {"dataSource": "data"}}, "phospho": {"pca": {"dataSource": "data"}}})
    assert all(options["data_source"] == "data" for name, _, options in figures if name == "qc_pca_plot")


def test_the_coverage_table_names_what_it_counts() -> None:
    import pandas as pd

    frame = pd.DataFrame({"Sample": ["WT_1"], "Number": [10]})
    assert "Identified proteins (n)" in report._table_to_html("Coverage", frame, "proteins")
    assert "Identified features (n)" in report._table_to_html("Coverage", frame)
