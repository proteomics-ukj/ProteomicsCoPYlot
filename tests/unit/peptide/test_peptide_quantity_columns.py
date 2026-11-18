"""Which peptide quantity column completeness reads.

Split out of test_small_bugs.py (2026-10-02)."""
from __future__ import annotations

import pandas as pd
import pytest

from app.services.peptide_quantity import resolve_peptide_quantity_column as _peptide_quantity_column


@pytest.mark.parametrize("header", ["Precursor.Quantity", "FG.Quantity"])
def test_imported_quantity_wins_over_group_measurement(header):
    from io import BytesIO
    from app.services.dataset_reader import read_dataframe

    report = f"{header}\tPG.Quantity\n0\t100\n"
    frame = read_dataframe("report.tsv", BytesIO(report.encode()))
    assert _peptide_quantity_column(frame, requested="PG.Quantity") == "Precursor.Quantity"
    assert frame["Precursor.Quantity"].tolist() == [0]


@pytest.mark.parametrize("header", ["PG.Quantity", "Genes.Quantity", "Intensity", "Quantity.Quality"])
def test_missing_precursor_quantity_is_reported(header):
    with pytest.raises(ValueError, match="FG.Quantity"):
        _peptide_quantity_column(pd.DataFrame({header: [100.0]}))

@pytest.mark.parametrize("spectronaut", [False, True], ids=["diann", "spectronaut"])
@pytest.mark.parametrize("level", ["peptide", "precursor"])
def test_coverage_and_missingness_use_imported_precursor_values(monkeypatch, spectronaut, level):
    from io import BytesIO
    from types import SimpleNamespace
    from app.services import completeness_tools, peptide_tools, plot_images
    from app.services.dataset_reader import read_dataframe

    headers = (["R.FileName", "PEP.StrippedSequence", "EG.PrecursorId", "FG.Quantity"] if spectronaut
               else ["Run", "Stripped.Sequence", "Precursor.Id", "Precursor.Quantity"])
    rows = ["r1\tAAA\tAAA2\t5\t100", "r1\tAAA\tAAA3\t0\t100",
            "r2\tAAA\tAAA2\t0\t100", "r2\tAAA\tAAA3\t0\t100"]
    raw = "\n".join(["\t".join([*headers, "PG.Quantity"]), *rows])
    frame = read_dataframe("report.tsv", BytesIO(raw.encode()))
    meta = pd.DataFrame({"sample": ["r1", "r2"], "condition": ["A", "B"]})
    monkeypatch.setattr(peptide_tools, "get_peptide_frame", lambda **kwargs: frame)
    monkeypatch.setattr(peptide_tools, "get_current_dataset", lambda _: None)
    monkeypatch.setattr(peptide_tools, "get_peptide_metadata", lambda: SimpleNamespace(frame=meta))
    coverage, _, _ = plot_images._peptide_coverage_summary(level=level)
    matrix = completeness_tools.completeness_peptide_matrix(level)
    assert coverage.Value.tolist() == [1, 0]
    assert matrix.notna().sum().tolist() == [1, 0]
    assert len(matrix) == (1 if level == "peptide" else 2)

def test_spectronaut_labeled_sequence_and_charge_identify_precursors():
    from app.services.dataset_reader import rename_cols
    from app.services.plot_images import peptide_precursor_column

    raw = pd.DataFrame({"FG.LabeledSequence": ["_AAA_", "_AAA_", "_AAA_", "_A[Ox]AA_", None, "_AAA_"],
                        "FG.Charge": [3, 2, 2, 2, 2, None]})
    frame = rename_cols(raw)
    ids = frame[peptide_precursor_column(frame)]
    assert ids.iloc[:4].tolist() == ["_AAA_/3", "_AAA_/2", "_AAA_/2", "_A[Ox]AA_/2"]
    assert ids.iloc[4:].isna().all()
    raw["EG.PrecursorId"] = "provided"
    direct = rename_cols(raw)
    assert direct[peptide_precursor_column(direct)].eq("provided").all()
