"""The coverage table per condition: mean, median and SD of the per-sample counts, and the same over all samples."""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import pytest

from proteomicscopylot import notebook_tables as pcp


def _data() -> tuple[pd.DataFrame, pd.DataFrame]:
    # Quantified features per sample: a1 4, a2 2, a3 3, b1 1, b2 1.
    data = pd.DataFrame({
        "ProteinNames": ["P1", "P2", "P3", "P4"],
        "a1": [1.0, 1.0, 1.0, 1.0], "a2": [1.0, 1.0, np.nan, np.nan], "a3": [1.0, 1.0, 1.0, np.nan],
        "b1": [1.0, np.nan, np.nan, np.nan], "b2": [np.nan, 1.0, np.nan, np.nan],
    })
    meta = pd.DataFrame({"sample": ["a1", "a2", "a3", "b1", "b2"], "condition": ["A", "A", "A", "B", "B"]})
    return data, meta


def test_each_condition_and_all_samples_get_mean_median_and_sd() -> None:
    table = pcp.qc_coverage_table(*_data(), summary=True)
    assert table.columns.tolist() == ["condition", "mean", "median", "sd", "min", "max"]
    assert table["condition"].tolist() == ["A", "B", "All samples"]
    assert table["mean"].tolist() == [3.0, 1.0, 2.2]
    assert table["median"].tolist() == [3.0, 1.0, 2.0]
    assert table["sd"].tolist() == [1.0, 0.0, 1.3]  # sample SD of (4, 2, 3), (1, 1) and all five
    assert table["min"].tolist() == [2, 1, 1] and table["max"].tolist() == [4, 1, 4]


def test_the_application_table_is_the_package_table() -> None:
    from app.services.table_functions import prepare_coverage_summary_df

    data, meta = _data()
    pd.testing.assert_frame_equal(prepare_coverage_summary_df(data, meta), pcp.qc_coverage_table(data, meta, summary=True))


def test_empty_filtered_matrix_keeps_zero_counts_for_all_samples() -> None:
    from app.services.table_functions import prepare_coverage_summary_df

    data, meta = _data()
    for table in (prepare_coverage_summary_df(data.iloc[:0], meta),
                  pcp.qc_coverage_table(data.iloc[:0], meta, summary=True)):
        assert table["condition"].tolist() == ["A", "B", "All samples"]
        assert (table[["mean", "median", "sd", "min", "max"]] == 0).all().all()


def test_workbook_preserves_formula_shaped_identifiers_as_text(monkeypatch) -> None:
    from openpyxl import load_workbook
    from app.api import plots_qc

    table = pd.DataFrame({"=1+1": ["=HYPERLINK(\"https://example.invalid\",\"protein\")"], "CV (%)": [2.5]})
    monkeypatch.setattr(plots_qc, "linear_qc_frame", lambda kind: pd.DataFrame())
    monkeypatch.setattr(plots_qc, "_get_sample_columns", lambda kind, frame: [])
    monkeypatch.setattr(plots_qc, "metadata_for_kind", lambda kind, frame, samples: pd.DataFrame())
    monkeypatch.setattr(plots_qc, "qc_feature_cvs", lambda frame, meta: table)
    response = plots_qc.qc_feature_cvs_route("protein", format="xlsx")
    sheet = load_workbook(io.BytesIO(response.body), data_only=False).active
    assert sheet["A1"].value == "=1+1" and sheet["A1"].data_type == "s"
    assert sheet["A2"].value == table.iloc[0, 0] and sheet["A2"].data_type == "s"
    assert sheet["B2"].value == 2.5


@pytest.mark.usefixtures("example_project")
def test_feature_cvs_can_be_downloaded_as_a_workbook() -> None:
    from app.api.plots_qc import qc_feature_cvs_route

    as_csv = qc_feature_cvs_route("protein")
    as_xlsx = qc_feature_cvs_route("protein", format="xlsx")
    assert as_csv.media_type == "text/csv"
    assert as_xlsx.media_type == "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
    assert as_xlsx.headers["content-disposition"].endswith('feature_cvs.xlsx"')
    expected = pd.read_csv(io.BytesIO(as_csv.body))
    pd.testing.assert_frame_equal(pd.read_excel(io.BytesIO(as_xlsx.body)), expected, check_dtype=False)
