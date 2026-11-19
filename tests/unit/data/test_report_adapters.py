"""Reports of other tools the reader reshapes (D34, 2026-11-13).

Small tables with the column layouts of the PRIDE examples the adapters were
built on: Spectronaut's peptide pivot (PXD071383), long report (PXD077951) and
PTM site report (PXD071048), MaxQuant's proteinGroups.txt and evidence.txt
(PXD085265), and FragPipe's combined protein, ion and peptide tables
(PXD076585; D37, 2026-11-18).
"""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import pytest

from app.services.dataset_reader import adapt_report, holds_protein_and_peptide_level, read_dataframe


def _read(frame: pd.DataFrame, *, peptide: bool = False, name: str = "report.tsv") -> pd.DataFrame:
    return read_dataframe(name, io.BytesIO(frame.to_csv(sep="\t", index=False).encode()), peptide=peptide)


# --- Spectronaut peptide pivot and long report ----------------------------------

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


def test_portable_notebook_reader_rejects_the_same_unsupported_inputs(tmp_path):
    from app.services.notebook_preprocessing import preprocessing_source

    frame, message = _site_pivot().assign(**{"PTM.ModificationTitle": "Oxidation (M)"}), "no phosphorylation sites"
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


# --- FragPipe ----------------------------------------------------------------------

def _fragpipe_ions() -> pd.DataFrame:
    return pd.DataFrame({
        "Peptide Sequence": ["PEPTIDEK", "MAGICR", "SAMPLER"],
        "Modified Sequence": ["PEPTIDEK", "n[42.0106]M[15.9949]AGIC[57.0215]R", "S[79.9663]AM[31.9898]PLER"],
        "Charge": [2, 3, 2],
        "Protein": ["sp|P1|ONE_HUMAN", "sp|P1|ONE_HUMAN", "contam_sp|P00761|TRYP_PIG"],
        "Protein ID": ["P1", "P1", "P00761"], "Entry Name": ["ONE_HUMAN", "ONE_HUMAN", "TRYP_PIG"], "Gene": ["G1", "G1", np.nan],
        "Mapped Proteins": [np.nan, "sp|P9|NINE_HUMAN", np.nan],
        "ctrl_1 Spectral Count": [2, 1, 0], "treat_1 Spectral Count": [1, 0, 3],
        "ctrl_1 Apex Retention Time": [600.0, 660.0, np.nan], "treat_1 Apex Retention Time": [606.0, np.nan, 720.0],
        "ctrl_1 Intensity": [1000.0, 2000.0, 0.0], "treat_1 Intensity": [1100.0, 0.0, 3100.0],
        "ctrl_1 Match Type": ["MS/MS", "MBR", "unmatched"], "treat_1 Match Type": ["MS/MS", "unmatched", "MS/MS"],
    })


def test_a_fragpipe_ion_table_becomes_a_long_report():
    long = _read(_fragpipe_ions(), peptide=True, name="combined_ion.tsv")
    assert long["File.Name"].tolist() == ["ctrl_1", "ctrl_1", "treat_1", "treat_1"]  # a quantity of 0 is no row
    assert long["Stripped.Sequence"].tolist() == ["PEPTIDEK", "MAGICR", "PEPTIDEK", "SAMPLER"]
    assert long["Precursor.Quantity"].tolist() == [1000, 2000, 1100, 3100]
    assert long["RT"].tolist() == [10.0, 11.0, 10.1, 12.0]  # seconds in the table, minutes here
    assert long["Match.Type"].tolist() == ["MS/MS", "MBR", "MS/MS", "MS/MS"]  # match between runs is kept and named
    # the masses the modification plot counts get their names, any other mass stays as written
    assert long["Modified.Sequence"].tolist()[1] == "n[Acetyl]M[Oxidation]AGIC[Carbamidomethyl]R"
    assert long["Modified.Sequence"].tolist()[3] == "S[Phospho]AM[31.9898]PLER"
    assert long["Precursor.Id"].tolist()[:2] == ["PEPTIDEK2", "n[Acetyl]M[Oxidation]AGIC[Carbamidomethyl]R3"]
    assert long["ProteinNames"].tolist() == ["ONE_HUMAN", "ONE_HUMAN", "ONE_HUMAN", "contam_sp|P00761|TRYP_PIG"]
    assert long["ProteinAccessions"].tolist()[0] == "P1" and long["GeneNames"].tolist()[0] == "G1"
    assert not [column for column in long.columns if "Spectral Count" in column or column.startswith(("ctrl_1 ", "treat_1 "))]


def test_the_peptide_table_uses_intensity_not_maxlfq_and_has_no_precursor():
    table = _fragpipe_ions().drop(columns=["Modified Sequence", "Charge", "ctrl_1 Apex Retention Time", "treat_1 Apex Retention Time"])
    table["ctrl_1 MaxLFQ Intensity"], table["treat_1 MaxLFQ Intensity"] = [7.0, 0.0, 0.0], [8.0, 0.0, 9.0]
    long = _read(table, peptide=True, name="combined_peptide.tsv")
    assert long["Precursor.Quantity"].tolist() == [1000, 2000, 1100, 3100]
    assert not {"Modified.Sequence", "Precursor.Id", "RT"} & set(long.columns)
    assert not [column for column in long.columns if "MaxLFQ" in column]


def test_the_fragpipe_columns_survive_the_detection_only_reading():
    from app.services.peptide_tools import _detection_column

    kept = [column for column in _fragpipe_ions().assign(**{"ctrl_1 MaxLFQ Intensity": 1.0}).columns if _detection_column(column)]
    assert kept == ["Peptide Sequence", "Modified Sequence", "Charge", "ctrl_1 Intensity", "treat_1 Intensity"]
    long = _read(_fragpipe_ions()[kept], peptide=True)
    assert {"File.Name", "Stripped.Sequence", "Precursor.Id", "Precursor.Quantity"}.issubset(long.columns) and len(long) == 4


def _fragpipe_proteins() -> pd.DataFrame:
    return pd.DataFrame({
        "Protein": ["sp|P1|ONE_HUMAN", "contam_sp|P00761|TRYP_PIG", "custom_entry", "contam_sp|P02769|ALBU_BOVIN", "contam_sp|P04264|K2C1_HUMAN"],
        "Protein ID": ["P1", "P00761", "custom_entry", "P02769", "P04264"],
        "Entry Name": ["ONE_HUMAN", "TRYP_PIG", np.nan, "ALBU_BOVIN", "K2C1_HUMAN"], "Gene": ["G1", np.nan, np.nan, "ALB", "KRT1"],
        "ctrl_1 Spectral Count": [3, 2, 1, 1, 1], "ctrl_1 Intensity": [10.0, 20.0, 30.0, 1.0, 1.0],
        "ctrl_1 MaxLFQ Intensity": [11.0, 0.0, 31.0, 1.0, 1.0],
        "Indistinguishable Proteins": ["sp|P2|TWO_HUMAN, tr|A0A9|A0A9_HUMAN, sp|P1-2|ONE_HUMAN", np.nan, np.nan, "sp|P02768|ALBU_HUMAN",
                                       "contam_sp|P35908|K22E_HUMAN"],
    })


def test_a_fragpipe_protein_table_is_named_by_its_group_and_keeps_the_contaminant_marker():
    from proteomicscopylot._core import upstream_contaminant_mask

    read = _read(_fragpipe_proteins(), name="combined_protein.tsv")
    assert read["ProteinNames"].tolist() == [
        "ONE_HUMAN;TWO_HUMAN;A0A9_HUMAN", "contam_sp|P00761|TRYP_PIG", "custom_entry",
        "contam_sp|P02769|ALBU_BOVIN;ALBU_HUMAN", "contam_sp|P04264|K2C1_HUMAN;contam_sp|P35908|K22E_HUMAN"]
    # an isoform shares its entry name (listed once) and keeps its own accession
    assert read["ProteinAccessions"].tolist() == ["P1;P2;A0A9;P1-2", "P00761", "custom_entry", "P02769;P02768", "P04264;P35908"]
    assert read["Protein"].tolist()[0] == "sp|P1|ONE_HUMAN" and read["Protein ID"].tolist()[0] == "P1"
    assert read["ctrl_1 MaxLFQ Intensity"].tolist() == [11.0, 0.0, 31.0, 1.0, 1.0] and len(read) == 5
    # a group is a contaminant only when every member is: bovine albumin grouped with the human one stays
    assert upstream_contaminant_mask(read).tolist() == [False, True, False, False, True]
    assert "ProteinNames" not in adapt_report(_fragpipe_proteins(), peptide=True).columns  # not a peptide table: left alone


def test_the_rt_plot_says_that_a_fragpipe_report_has_no_predicted_retention_time(api):
    client, _ = api
    _upload(client, "/api/datasets/peptide-file", _fragpipe_ions().assign(**{"ctrl_1 Match Type": ["MS/MS", "MBC", "unmatched"]}),
            "combined_ion.tsv")
    response = client.get("/api/plots/peptide/rt.png")
    assert response.status_code == 400 and "no predicted retention time" in response.json()["detail"]
    from app.services.peptide_tools import get_peptide_frame
    assert get_peptide_frame()["Match.Type"].tolist() == ["MS/MS", "MBC", "MS/MS", "MS/MS"]


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
    (["KD_Rad_1 MaxLFQ Intensity", "KD_1 MaxLFQ Intensity"], ["KD_Rad_1", "KD_1"]),  # FragPipe
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


def test_a_long_report_loaded_as_the_protein_table_fills_the_peptide_report(api):
    client, loaded = api
    preview = _upload(client, "/api/datasets/upload", _spectronaut_long(), "long.tsv", kind="protein")
    assert preview["rows"] == 2 and loaded("peptide") == "long.tsv"
    overview = client.get("/api/peptide/overview").json()
    assert overview["rows"] == 6 and overview["availableSamples"] == ["ctrl_1", "treat_1"]


def test_loaded_as_the_peptide_report_it_fills_the_protein_table(api):
    client, loaded = api
    _upload(client, "/api/datasets/peptide-file", _spectronaut_long(), "long.tsv")
    assert loaded("protein") == "long.tsv" and loaded("peptide") == "long.tsv"
    current = client.get("/api/datasets/current").json()
    assert current["protein"]["rows"] == 2


def test_another_protein_table_is_not_replaced(api):
    client, loaded = api
    _upload(client, "/api/datasets/upload", _protein_groups(), "proteinGroups.txt", kind="protein")
    _upload(client, "/api/datasets/peptide-file", _spectronaut_long(), "long.tsv")
    assert loaded("protein") == "proteinGroups.txt" and loaded("peptide") == "long.tsv"


def test_a_site_report_fills_nothing_else(api):
    client, loaded = api
    _upload(client, "/api/datasets/upload", _site_pivot(), "sites.tsv", kind="phospho")
    assert loaded("phospho") == "sites.tsv" and loaded("peptide") is None and loaded("protein") is None


def test_loading_it_again_as_the_other_level_keeps_what_is_there(api):
    client, loaded = api
    from app.services.dataset_store import get_current_dataset
    _upload(client, "/api/datasets/upload", _spectronaut_long(), "long.tsv", kind="protein")
    protein = get_current_dataset("protein")
    _upload(client, "/api/datasets/peptide-file", _spectronaut_long(), "long.tsv")
    assert get_current_dataset("protein") is protein and loaded("peptide") == "long.tsv"
