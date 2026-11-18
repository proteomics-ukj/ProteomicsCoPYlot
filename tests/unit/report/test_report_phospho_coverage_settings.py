"""QC counts and localization coverage use their own settings in reports."""
from __future__ import annotations

import importlib
from types import SimpleNamespace

import pandas as pd

from app.schemas.summary import SummaryReportContext, SummaryReportRequest

report = importlib.import_module("app.services.report_function")


def _render(monkeypatch, qc_settings):
    annotation = SimpleNamespace(
        metadata=pd.DataFrame({"sample": ["A", "B"], "condition": ["A", "B"]}),
        source_data=pd.DataFrame({"value": [1, 2, 3]}), filtered_data=pd.DataFrame({"value": [1, 2]}),
        filter_config=SimpleNamespace(method="percentage", mode="overall", minPercentage=75, normalization="median"),
        imputation=SimpleNamespace(mode="numerical", downShift=1.8, width=0.3, seed=69, sampleWise=True),
    )
    monkeypatch.setattr(report, "get_annotation", lambda kind: annotation if kind == "phospho" else None)
    monkeypatch.setattr(report, "get_current_dataset", lambda kind: object() if kind == "phospho" else None)
    calls = []

    def renderer(name):
        def render(*args, **kwargs):
            calls.append((name, kwargs))
            if name.endswith("_table"):
                return pd.DataFrame({"Label": ["x"], "Value": [1]})
            return "<div></div>" if name.endswith("_html") else b"png"
        return render

    for name in [n for n in vars(report) if n.startswith(("qc_", "peptide_")) or n in (
            "completeness_missing_value_plot", "phospho_coverage_png", "phospho_sty_png", "phosphosite_plot_png", "volcano_html")]:
        if callable(getattr(report, name)):
            monkeypatch.setattr(report, name, renderer(name))
    context = SummaryReportContext(qcSettings={"phospho": qc_settings})
    result = report.report_function(SummaryReportRequest(title="t", author="a", reportContext=context))
    return result.html, calls


def test_coverage_plots_use_their_own_page_settings(monkeypatch):
    html, calls = _render(monkeypatch, {
        "coverage": {"randomSeed": 5, "includeId": False, "summary": True},
        "phosphoCoverage": {"includeId": True, "conditions": ["B"], "colorClassI": "#111111",
                            "colorNotClassI": "#222222", "randomSeed": 7, "mode": "Summary"},
    })
    qc = [kw for name, kw in calls if name == "qc_coverage_plot"]
    localized = [kw for name, kw in calls if name == "phospho_coverage_png"]
    assert len(qc) == len(localized) == 1
    assert qc[0]["kind"] == "phospho" and qc[0]["summary"] is True
    assert qc[0]["random_seed"] == 5 and qc[0]["include_id"] is False
    call = localized[0]
    assert call["mode"] == "Summary"
    assert call["include_id"] is True and call["conditions"] is None
    assert (call["color_class_i"], call["color_not_class_i"]) == ("#111111", "#222222")
    assert call["random_seed"] == 7
    assert html.count(report.CLASS_I_SENTENCE) == 1
    assert "Coverage Plot (per condition)" in html
    assert "Coverage Plot (localization probability)" in html
    assert "seed" not in html.lower()


def test_page_defaults_apply_without_saved_settings(monkeypatch):
    _, calls = _render(monkeypatch, {})
    localized = [kw for name, kw in calls if name == "phospho_coverage_png"]
    assert len(localized) == 1
    assert localized[0]["mode"] == "Normal"
    assert localized[0]["include_id"] is False and localized[0]["conditions"] is None
    assert (localized[0]["color_class_i"], localized[0]["color_not_class_i"]) == ("", "")


def test_sty_uses_saved_size_and_has_its_own_number(monkeypatch):
    html, calls = _render(monkeypatch, {"sty": {"widthCm": 18, "heightCm": 12, "dpi": 120}})
    settings = next(kw for name, kw in calls if name == "phospho_sty_png")
    assert (settings["width_cm"], settings["height_cm"], settings["dpi"]) == (18, 12, 120)
    assert "<h2 id='s-1-4'>1.4 STY Plot</h2>" in html
    assert "<h2 id='s-1-5'>1.5 Missing Value Plot</h2>" in html
