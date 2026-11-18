"""Reports of other tools the reader reshapes (D34, 2026-11-13).

Small tables with the column layouts of the PRIDE examples the adapters were
built on: a Skyline transition export (PXD053085), Spectronaut's peptide pivot
(PXD071383), long report (PXD077951) and PTM site report (PXD071048), and
MaxQuant's proteinGroups.txt and evidence.txt (PXD085265).
"""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import pytest

from app.services.dataset_reader import adapt_report, holds_protein_and_peptide_level, read_dataframe


def _read(frame: pd.DataFrame, *, peptide: bool = False, name: str = "report.tsv") -> pd.DataFrame:
    return read_dataframe(name, io.BytesIO(frame.to_csv(sep="\t", index=False).encode()), peptide=peptide)


# --- Skyline -------------------------------------------------------------------

def _skyline() -> pd.DataFrame:
    rows = []
    for protein, peptide, modified, charge in [
        ("sp|P11111|ONE_HUMAN", "PEPTIDEK", "PEPTIDEK", 2),
        ("sp|P11111|ONE_HUMAN", "MAGICR", "M[+16]AGICR", 3),
        ("sp|P22222|TWO_HUMAN", "SAMPLER", "SAMPLER", 2),
    ]:
        for run, scale in [("a_1.wiff", 1.0), ("b_1.wiff", 2.0)]:
            for fragment, area in [("y3", 100.0), ("y4", 50.0)]:
                rows.append({"Protein Name": protein, "Peptide Sequence": peptide, "Peptide Modified Sequence": modified,
                             "Precursor Charge": charge, "Precursor Mz": 500.5, "Fragment Ion": fragment, "Product Charge": 1,
                             "Product Mz": 400.2, "Isotope Label Type": "light", "File Name": run, "Area": area * scale,
                             "Truncated": False, "Detection Q Value": 0.01})
    frame = pd.DataFrame(rows)
    # One protein is not measured in the second run, and a heavy standard is not a measurement of the sample.
    frame.loc[(frame["Protein Name"] == "sp|P22222|TWO_HUMAN") & (frame["File Name"] == "b_1.wiff"), "Area"] = 0
    heavy = frame.iloc[[0]].assign(**{"Isotope Label Type": "heavy", "Area": 9999.0})
    return pd.concat([frame, heavy], ignore_index=True)


def test_a_skyline_export_is_the_protein_matrix_when_loaded_as_a_table():
    matrix = _read(_skyline(), name="skyline.tsv")
    assert matrix.columns.tolist() == ["ProteinNames", "ProteinAccessions", "a_1.wiff", "b_1.wiff"]
    one, two = matrix.set_index("ProteinNames").loc[["ONE_HUMAN", "TWO_HUMAN"]].to_dict("records")
    assert (one["a_1.wiff"], one["b_1.wiff"]) == (300.0, 600.0)  # both precursors, all light transitions
    assert two["a_1.wiff"] == 150.0 and np.isnan(two["b_1.wiff"])  # zero areas are not a measurement
    assert matrix["ProteinAccessions"].tolist() == ["P11111", "P22222"]


def test_the_same_export_is_the_peptide_report_when_loaded_as_one():
    long = _read(_skyline(), peptide=True, name="skyline.tsv")
    assert {"File.Name", "ProteinNames", "Stripped.Sequence", "Modified.Sequence", "Precursor.Id", "Precursor.Quantity"} <= set(long.columns)
    assert len(long) == 6  # three precursors in two runs
    oxidised = long[long["Modified.Sequence"] == "M[+16]AGICR"]
    assert oxidised["Stripped.Sequence"].unique().tolist() == ["MAGICR"]
    assert oxidised["Precursor.Id"].unique().tolist() == ["M[+16]AGICR3"]
    first = long[(long["Modified.Sequence"] == "PEPTIDEK") & (long["File.Name"] == "a_1.wiff")]
    assert first["Precursor.Quantity"].tolist() == [150.0]  # the heavy row is left out


def _skyline_with_undetected_peaks() -> pd.DataFrame:
    frame = _skyline()
    one_b = (frame["Protein Name"] == "sp|P11111|ONE_HUMAN") & (frame["File Name"] == "b_1.wiff") & (frame["Isotope Label Type"] == "light")
    # In run b the first precursor is not detected (q above 0.01) and one transition of the second is truncated.
    frame.loc[one_b & (frame["Peptide Sequence"] == "PEPTIDEK"), "Detection Q Value"] = 0.2
    frame.loc[one_b & (frame["Peptide Sequence"] == "MAGICR") & (frame["Fragment Ion"] == "y4"), "Truncated"] = True
    # A peak without a q-value is kept; decoy proteins and iRT standards are not sample proteins.
    frame.loc[(frame["Protein Name"] == "sp|P22222|TWO_HUMAN") & (frame["File Name"] == "a_1.wiff"), "Detection Q Value"] = np.nan
    frame["Standard Type"] = ""
    decoy = frame.iloc[[1]].assign(**{"Protein Name": "Decoys"})
    irt = frame.iloc[[1]].assign(**{"Protein Name": "Biognosys|iRT-Kit", "Standard Type": "iRT"})
    # A protein whose peaks were detected in no run is not listed.
    absent = frame.iloc[[0, 1, 4, 5]].assign(**{"Protein Name": "sp|P33333|THREE_HUMAN", "Detection Q Value": 0.9})
    return pd.concat([frame, decoy, irt, absent], ignore_index=True)


def test_skyline_peaks_that_were_not_detected_are_not_measurements():
    """The defaults of MSstats' SkylinetoMSstatsFormat: q above 0.01, truncated peaks, decoys and iRT are left out."""
    matrix = _read(_skyline_with_undetected_peaks(), name="skyline.tsv")
    assert matrix["ProteinNames"].tolist() == ["ONE_HUMAN", "TWO_HUMAN"]
    one, two = matrix.set_index("ProteinNames").loc[["ONE_HUMAN", "TWO_HUMAN"]].to_dict("records")
    assert one["a_1.wiff"] == 300.0
    assert one["b_1.wiff"] == 200.0  # only the y3 transition of the second precursor: 100 * 2
    assert two["a_1.wiff"] == 150.0  # no q-value: kept
    notes = matrix.attrs["loading_notes"]
    assert any(note.startswith("2 transition rows of decoy proteins or iRT standards") for note in notes)
    assert any(note.startswith("1 of 14 measured transition areas are truncated") for note in notes)
    assert any(note.startswith("6 of 14 measured transition areas have a Detection Q Value above 0.01") for note in notes)
    assert notes[-1] == "1 proteins without a detected peak in any run are not listed."

    long = _read(_skyline_with_undetected_peaks(), peptide=True, name="skyline.tsv")
    undetected = long[(long["ProteinNames"] == "ONE_HUMAN") & (long["Modified.Sequence"] == "PEPTIDEK") & (long["File.Name"] == "b_1.wiff")]
    assert undetected["Precursor.Quantity"].isna().all()
    assert long.attrs["loading_notes"] == notes[:-1]


def test_the_annotation_states_what_a_skyline_import_left_out():
    from app.services.annotation_processor import _table_shape_warnings
    warnings = _table_shape_warnings(_read(_skyline_with_undetected_peaks(), name="skyline.tsv"))
    assert any("Detection Q Value above 0.01" in warning for warning in warnings)


# --- Spectronaut peptide pivot and long report ----------------------------------

@pytest.mark.parametrize("peptide", [False, True])
@pytest.mark.parametrize("label", ["heavy", "", None])
def test_skyline_without_light_measurements_is_rejected(peptide, label):
    frame = _skyline().assign(**{"Isotope Label Type": label})
    with pytest.raises(ValueError, match="no light-labelled measurements"):
        _read(frame, peptide=peptide)


def test_skyline_without_a_label_column_still_loads():
    frame = _skyline()
    frame = frame.loc[frame["Isotope Label Type"] == "light"].drop(columns="Isotope Label Type")
    assert _read(frame).set_index("ProteinNames").loc["ONE_HUMAN", "a_1.wiff"] == 300.0


def _spectronaut_pivot() -> pd.DataFrame:
    return pd.DataFrame({
        "PG.ProteinGroups": ["P1", "P1", "P2"], "PG.Genes": ["G1", "G1", "G2"],
        "EG.PrecursorId": ["_PEPTIDEK_.2", "_M[Oxidation (M)]AGICR_.3", "_SAMPLER_.2"],
        "[1] ctrl 1.raw.EG.TotalQuantity (Settings)": [10.0, 20.0, "Filtered"],
        "[2] treat 1.raw.EG.TotalQuantity (Settings)": [11.0, np.nan, 31.0],
        "[1] ctrl 1.raw.EG.ApexRT": [5.0, 6.0, np.nan], "[2] treat 1.raw.EG.ApexRT": [5.1, np.nan, 7.1],
        "[1] ctrl 1.raw.EG.RTPredicted": [5.2, 6.2, np.nan], "[2] treat 1.raw.EG.RTPredicted": [5.3, np.nan, 7.3],
        "[1] ctrl 1.raw.EG.Qvalue": [0.001, 0.002, 1.0], "[2] treat 1.raw.EG.Qvalue": [0.001, 1.0, 0.003],
    })


def test_a_spectronaut_peptide_pivot_becomes_a_long_report():
    long = _read(_spectronaut_pivot(), peptide=True)
    assert long["File.Name"].tolist() == ["[1] ctrl 1.raw", "[1] ctrl 1.raw", "[2] treat 1.raw", "[2] treat 1.raw"]
    assert long["Precursor.Quantity"].tolist() == [10.0, 20.0, 11.0, 31.0]  # no row where nothing was quantified
    assert long["Stripped.Sequence"].tolist() == ["PEPTIDEK", "MAGICR", "PEPTIDEK", "SAMPLER"]
    assert long["Modified.Sequence"].iloc[1] == "_M[Oxidation (M)]AGICR_"
    assert long["RT"].tolist() == [5.0, 6.0, 5.1, 7.1] and long["Predicted.RT"].tolist() == [5.2, 6.2, 5.3, 7.3]
    assert not [column for column in long.columns if "Qvalue" in column or "TotalQuantity" in column]


def test_the_same_pivot_stays_a_matrix_as_a_table():
    table = _read(_spectronaut_pivot())
    assert len(table) == 3 and "[1] ctrl 1.raw.EG.TotalQuantity (Settings)" in table.columns


def _spectronaut_long() -> pd.DataFrame:
    rows = []
    for run, condition, scale in [("ctrl_1", "ctrl", 1.0), ("treat_1", "treat", 2.0)]:
        for group, quantity, sequence, modified in [("P1", 100.0, "PEPTIDEK", "_PEPTIDEK_"), ("P1", 100.0, "MAGICR", "_MAGICR_"),
                                                    ("P2;P3", 40.0, "SAMPLER", "_SAMPLER_")]:
            rows.append({"R.Condition": condition, "R.FileName": run, "PG.ProteinGroups": group, "PG.ProteinAccessions": group,
                         "PG.Qvalue": 0.001, "PG.Quantity": quantity * scale, "PEP.StrippedSequence": sequence,
                         "EG.ModifiedSequence": modified, "EG.ApexRT": 10.0, "EG.RTPredicted": 10.1,
                         "EG.TotalQuantity (Settings)": 7.0 * scale, "FG.Charge": 2})
    return pd.DataFrame(rows)


def test_a_spectronaut_long_report_is_the_protein_matrix_as_a_table():
    matrix = _read(_spectronaut_long())
    assert len(matrix) == 2
    assert matrix["ctrl_1.PG.Quantity"].tolist() == [100.0, 40.0] and matrix["treat_1.PG.Quantity"].tolist() == [200.0, 80.0]
    assert "ProteinNames" in matrix.columns and "PG.Qvalue" not in matrix.columns


def test_and_the_peptide_report_as_one():
    long = _read(_spectronaut_long(), peptide=True)
    assert len(long) == 6 and long["Precursor.Quantity"].tolist()[:3] == [7.0, 7.0, 7.0]
    assert {"File.Name", "Stripped.Sequence", "Modified.Sequence", "RT", "Predicted.RT"} <= set(long.columns)


def test_which_reports_hold_both_levels():
    assert holds_protein_and_peptide_level(_skyline().columns)
    assert holds_protein_and_peptide_level(_spectronaut_long().columns)
    assert not holds_protein_and_peptide_level(_spectronaut_pivot().columns)
    assert not holds_protein_and_peptide_level(["Protein.Group", "Run", "Precursor.Id", "Precursor.Quantity"])


# --- Spectronaut PTM site report -------------------------------------------------

def _site_pivot(**extra) -> pd.DataFrame:
    return pd.DataFrame({
        "PTM.ProteinId": ["P1", "P1", "P2", "P2"], "PTM.CollapseKey": ["P1_S30_M1", "P1_T35_M2", "P2_Y7_M1", "P2_M1_M1"],
        "PTM.Multiplicity": [1, 2, 1, 1], "PTM.ModificationTitle": ["Phospho", "Phospho (STY)", "Phospho", "Oxidation (M)"],
        "PTM.SiteAA": ["S", "T", "Y", "M"], "PTM.SiteLocation": [30, 35, 7, 1],
        "PTM.FlankingRegion": ["YRPPPLDSIRLPRYV", "LDSIRLPTYVAAAAA", "______MYAAAAAAA", "_______MAAAAAAA"],
        "[1] a.raw.PTM.Quantity": [1.0, 2.0, 3.0, 4.0], "[2] b.raw.PTM.Quantity": [1.5, np.nan, 3.5, 4.5],
        "[1] a.raw.PTM.Stoichiometry": [np.nan] * 4, "[2] b.raw.PTM.Stoichiometry": [np.nan] * 4, **extra})


def test_spectronauts_site_report_gets_the_apps_phosphosite_columns():
    sites = _read(_site_pivot())
    assert sites["PTM_Collapse_key"].tolist() == ["P1~P1_S30_M1", "P1~P1_T35_M2", "P2~P2_Y7_M1"]  # phospho rows only
    assert sites["Protein_group"].tolist() == ["P1", "P1", "P2"] and sites["Gene_group"].tolist() == ["P1", "P1", "P2"]
    assert sites["UPD_seq"].tolist() == ["YRPPPLDs*IRLPRYV", "LDSIRLPt*YVAAAAA", "______My*AAAAAAA"]
    assert sites["PTM_localization"].isna().all()  # exported without a site probability
    assert [column for column in sites.columns if column.endswith("PTM.Quantity")] == ["[1] a.raw.PTM.Quantity", "[2] b.raw.PTM.Quantity"]
    assert not [column for column in sites.columns if "Stoichiometry" in column]


@pytest.mark.parametrize("title", ["Oxidation (M)", "", None])
@pytest.mark.parametrize("long_report", [False, True])
def test_spectronaut_without_phosphorylation_sites_is_rejected(title, long_report):
    frame = _site_pivot().assign(**{"PTM.ModificationTitle": title})
    if long_report:
        frame = frame.rename(columns={"[1] a.raw.PTM.Quantity": "PTM.Quantity"})
        frame["R.FileName"] = "run_a"
    with pytest.raises(ValueError, match="no phosphorylation sites"):
        _read(frame)


@pytest.mark.parametrize("report_type", ["skyline", "spectronaut"])
def test_portable_notebook_reader_rejects_the_same_unsupported_inputs(tmp_path, report_type):
    from app.services.notebook_preprocessing import preprocessing_source

    frame, message = (
        (_skyline().assign(**{"Isotope Label Type": "heavy"}), "no light-labelled measurements")
        if report_type == "skyline" else
        (_site_pivot().assign(**{"PTM.ModificationTitle": "Oxidation (M)"}), "no phosphorylation sites")
    )
    path = tmp_path / "report.tsv"
    frame.to_csv(path, sep="\t", index=False)
    namespace = {}
    exec(preprocessing_source().code, namespace)
    with pytest.raises(ValueError, match=message):
        namespace["read_original_table"](path)


def test_genes_and_site_probabilities_are_used_when_the_report_has_them():
    sites = _read(_site_pivot(**{"PG.Genes": ["AKT1;AKT2", "AKT1;AKT2", "MY_G", "MY_G"],
                                 "[1] a.raw.PTM.SiteProbability": [0.9, 0.4, 0.99, 1.0],
                                 "[2] b.raw.PTM.SiteProbability": [0.7, np.nan, 1.0, 1.0]}))
    assert sites["PTM_Collapse_key"].tolist() == ["P1~AKT1_S30_M1", "P1~AKT1_T35_M2", "P2~MY#G_Y7_M1"]
    assert sites["PTM_localization"].tolist() == [0.9, 0.4, 1.0]  # the site's highest probability


def test_a_long_site_report_is_pivoted_to_one_column_per_run():
    rows = []
    for run, scale in [("run_a", 1.0), ("run_b", 2.0)]:
        for key, residue, position, probability in [("P1_S30_M1", "S", 30, 0.8), ("P2_Y7_M1", "Y", 7, 0.6)]:
            rows.append({"R.FileName": run, "R.Condition": "x", "PG.Genes": "G" + key[1], "PTM.ProteinId": key[:2],
                         "PTM.CollapseKey": key, "PTM.ModificationTitle": "Phospho (STY)", "PTM.Multiplicity": 1,
                         "PTM.SiteLocation": position, "PTM.SiteAA": residue, "PTM.SiteProbability": probability * scale / 2,
                         "PTM.FlankingRegion": "AAAAAAA" + residue + "AAAAAAA", "PTM.Quantity": 10.0 * scale})
    sites = _read(pd.DataFrame(rows))
    assert sites["PTM_Collapse_key"].tolist() == ["P1~G1_S30_M1", "P2~G2_Y7_M1"]
    assert sites["run_a.PTM.Quantity"].tolist() == [10.0, 10.0] and sites["run_b.PTM.Quantity"].tolist() == [20.0, 20.0]
    assert sites["PTM_localization"].tolist() == [0.8, 0.6]


def test_the_collapsed_table_satisfies_the_phosphoprotein_collapse():
    from app.services.phosprot_tools import _REQUIRED_PHOSPHO_COLUMNS
    assert set(_REQUIRED_PHOSPHO_COLUMNS) <= set(_read(_site_pivot()).columns)


# --- MaxQuant --------------------------------------------------------------------

def _protein_groups() -> pd.DataFrame:
    return pd.DataFrame({
        "Protein IDs": ["P1", "REV__P2", "CON__P3", "P4"], "Majority protein IDs": ["P1", "REV__P2", "CON__P3", "P4"],
        "Gene names": ["G1", "", "", "G4"], "LFQ intensity a": [100, 0, 300, 400], "LFQ intensity b": [110, 0, 310, 0],
        "Only identified by site": ["", "", "", "+"], "Reverse": ["", "+", "", ""], "Potential contaminant": ["", "", "+", ""]})


def test_maxquant_decoys_and_site_only_groups_are_removed_on_loading_and_contaminants_stay():
    table = _read(_protein_groups(), name="proteinGroups.txt")
    assert table["ProteinNames"].tolist() == ["P1", "CON__P3"]
    assert table.attrs["decoy_rows_removed"] == 1
    assert table.attrs["loading_notes"][0].startswith("1 protein groups identified only by a modification site")


def test_the_annotation_says_how_many_rows_were_removed_on_loading():
    from app.services.annotation_processor import _table_shape_warnings
    warnings = _table_shape_warnings(_read(_protein_groups(), name="proteinGroups.txt"))
    assert any(warning.startswith("1 decoy rows") for warning in warnings)
    assert any(warning.startswith("1 protein groups identified only by a modification site") for warning in warnings)
    assert not any("Only identified by site)" in warning for warning in warnings)
    assert any(warning.startswith("1 rows look like decoy or contaminant entries") for warning in warnings)


def test_evidence_txt_reads_as_a_peptide_report():
    evidence = pd.DataFrame({
        "Sequence": ["PEPTIDEK", "MAGICR", "DECOYK"], "Modified sequence": ["_PEPTIDEK_", "_M(Oxidation (M))AGICR_", "_DECOYK_"],
        "Leading razor protein": ["P1", "P1", "REV__P9"], "Gene names": ["G1", "G1", ""], "Raw file": ["raw_01", "raw_02", "raw_01"],
        "Experiment": ["ctrl", "treat", "ctrl"], "Charge": [2, 3, 2], "Retention time": [10.0, 11.0, 12.0],
        "Intensity": [1000, 2000, 3000], "Reverse": ["", "", "+"]})
    long = _read(evidence, peptide=True, name="evidence.txt")
    assert len(long) == 2  # the decoy is gone
    assert long["File.Name"].tolist() == ["ctrl", "treat"]  # the experiment, as proteinGroups.txt names it
    assert long["Stripped.Sequence"].tolist() == ["PEPTIDEK", "MAGICR"]
    assert long["Precursor.Id"].tolist() == ["_PEPTIDEK_2", "_M(Oxidation (M))AGICR_3"]
    assert long["Precursor.Quantity"].tolist() == [1000, 2000] and long["RT"].tolist() == [10.0, 11.0]
    assert long["ProteinNames"].tolist() == ["P1", "P1"]


# --- what is left alone -----------------------------------------------------------

@pytest.mark.parametrize("peptide", [False, True])
def test_tables_without_an_adapter_pass_through(peptide):
    matrix = pd.DataFrame({"Protein.Group": ["P1"], "Genes": ["G1"], "run_a.raw": [1.0], "run_b.raw": [2.0]})
    report = pd.DataFrame({"File.Name": ["a", "b"], "Protein.Group": ["P1", "P1"], "Precursor.Id": ["PEPK2", "PEPK2"],
                           "Stripped.Sequence": ["PEPK", "PEPK"], "Precursor.Quantity": [1.0, 2.0]})
    for frame in (matrix, report):
        pd.testing.assert_frame_equal(adapt_report(frame, peptide), frame)


# --- protein and peptide runs of one search ---------------------------------------

@pytest.mark.parametrize("samples, runs", [
    (["[1] ctrl 1.raw.PG.Quantity", "[2] treat 1.raw.PG.Quantity"], ["[1] ctrl 1.raw", "[2] treat 1.raw"]),
    (["LFQ intensity ctrl_1", "LFQ intensity treat_1"], ["ctrl_1", "treat_1"]),
    (["ctrl_1.PG.Quantity", "treat_1.PG.Quantity"], ["ctrl_1", "treat_1"]),
])
def test_peptide_runs_are_matched_to_the_protein_samples_of_the_same_report(samples, runs):
    from app.services.peptide_tools import _sample_rename_map
    meta = pd.DataFrame({"sample": samples, "new_sample": ["ctrl_1", "treat_1"]})
    assert _sample_rename_map(runs, meta, "File.Name") == dict(zip(runs, ["ctrl_1", "treat_1"]))


def test_a_stem_two_samples_share_is_not_used():
    from app.services.peptide_tools import _sample_rename_map
    meta = pd.DataFrame({"sample": ["LFQ intensity x_7", "Intensity x_7"], "new_sample": ["a_1", "b_1"]})
    with pytest.raises(ValueError, match="ambiguous"):
        _sample_rename_map(["x_7"], meta, "File.Name")


# --- one upload, both levels --------------------------------------------------------

@pytest.fixture()
def api():
    from fastapi.testclient import TestClient

    from app.main import app
    from app.services.dataset_store import get_current_dataset, _package_state

    state = _package_state(None)
    before = dict(state)
    for kind in ("protein", "phospho", "phosprot", "peptide"):
        state[kind] = None
    with TestClient(app, base_url="http://localhost") as client:
        yield client, lambda kind: getattr(get_current_dataset(kind), "filename", None)
    state.clear()
    state.update(before)


def _upload(client, url: str, frame: pd.DataFrame, name: str, **form):
    body = frame.to_csv(sep="\t", index=False).encode()
    response = client.post(url, files={"file": (name, body, "text/tab-separated-values")}, data=form)
    assert response.status_code == 200, response.text
    return response.json()


def test_a_skyline_export_loaded_as_the_protein_table_fills_the_peptide_report(api):
    client, loaded = api
    preview = _upload(client, "/api/datasets/upload", _skyline(), "skyline.tsv", kind="protein")
    assert preview["rows"] == 2 and loaded("peptide") == "skyline.tsv"
    overview = client.get("/api/peptide/overview").json()
    assert overview["rows"] == 6 and overview["availableSamples"] == ["a_1.wiff", "b_1.wiff"]


def test_loaded_as_the_peptide_report_it_fills_the_protein_table(api):
    client, loaded = api
    _upload(client, "/api/datasets/peptide-file", _spectronaut_long(), "long.tsv")
    assert loaded("protein") == "long.tsv" and loaded("peptide") == "long.tsv"
    current = client.get("/api/datasets/current").json()
    assert current["protein"]["rows"] == 2


def test_another_protein_table_is_not_replaced(api):
    client, loaded = api
    _upload(client, "/api/datasets/upload", _protein_groups(), "proteinGroups.txt", kind="protein")
    _upload(client, "/api/datasets/peptide-file", _skyline(), "skyline.tsv")
    assert loaded("protein") == "proteinGroups.txt" and loaded("peptide") == "skyline.tsv"


def test_a_site_report_fills_nothing_else(api):
    client, loaded = api
    _upload(client, "/api/datasets/upload", _site_pivot(), "sites.tsv", kind="phospho")
    assert loaded("phospho") == "sites.tsv" and loaded("peptide") is None and loaded("protein") is None


def test_loading_it_again_as_the_other_level_keeps_what_is_there(api):
    client, loaded = api
    from app.services.dataset_store import get_current_dataset
    _upload(client, "/api/datasets/upload", _skyline(), "skyline.tsv", kind="protein")
    protein = get_current_dataset("protein")
    _upload(client, "/api/datasets/peptide-file", _skyline(), "skyline.tsv")
    assert get_current_dataset("protein") is protein and loaded("peptide") == "skyline.tsv"
