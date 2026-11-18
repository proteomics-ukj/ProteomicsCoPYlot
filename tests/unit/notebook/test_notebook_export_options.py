"""What a notebook export ships, and where its input comes from.

Two options, both off by default, decide it: data matrices of deterministic
steps, and the full imputed matrix. Without them an export holds the code, the
settings, the outputs and the values imputation drew; the input is rebuilt from
the original file, whose path is set in the first cell. With a matrix shipped
the notebook runs as it is and switches to the original file once its path is
filled in. The file is checked by content fingerprint (a mismatch warns), and
the last cell always compares the outputs byte for byte.
"""
from __future__ import annotations

import io
import json
import subprocess
import sys
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

from app.api.routes_annotations import generate_annotation
from app.schemas.annotation import AnnotationGenerateRequest
from app.schemas.data_tools import ImputationRunRequest
from app.services import data_tools
from app.services.dataset_reader import read_dataframe
from app.services.dataset_store import save_table_dataset
from app.services.notebook_archive import build_notebook_archive
from app.services.notebook_snapshots import snapshots
from app.services.session_context import get_session_id

SAMPLES = ["A1", "A2", "A3", "B1", "B2", "B3"]
CONDITIONS = [{"name": group, "columns": [f"{group}{i}" for i in range(1, 4)]} for group in "AB"]
RUNNER = """
import json, glob, sys, warnings
import IPython.display as D
D.display = lambda *a, **k: None
notebook = json.load(open(glob.glob('*.ipynb')[0]))
scope = {'__name__': '__main__'}
for cell in notebook['cells']:
    if cell['cell_type'] == 'code':
        source = ''.join(cell['source'])
        for old, new in json.loads(sys.argv[1]).items():
            source = source.replace(old, new)
        exec(compile(source, 'cell', 'exec'), scope)
"""


@pytest.fixture(autouse=True)
def empty_session():
    from app.services.annotation_store import reset_all_annotations
    from app.services.dataset_store import reset_all_datasets
    reset_all_annotations()
    reset_all_datasets()
    yield
    reset_all_annotations()
    reset_all_datasets()


@pytest.fixture
def analysis(tmp_path):
    values = np.random.default_rng(5).uniform(8, 16, (40, 6))
    values.flat[np.random.default_rng(6).choice(values.size, 30, replace=False)] = np.nan
    raw = pd.DataFrame(values, columns=SAMPLES)
    raw["ProteinNames"] = [f"P{i}" for i in range(40)]
    path = tmp_path / "protein.tsv"
    raw.to_csv(path, sep="\t", index=False)
    with path.open("rb") as stream:
        save_table_dataset(path.name, "protein", read_dataframe(path.name, stream), source_path=str(path))
    generate_annotation(AnnotationGenerateRequest(kind="protein", isLog2Transformed=True, conditions=CONDITIONS))
    data_tools.run_imputation(ImputationRunRequest(kind="protein"))
    return path


def _export(token, directory, **options):
    name, content = build_notebook_archive(snapshots.get(get_session_id(), token), **options)
    directory.mkdir()
    ZipFile(io.BytesIO(content)).extractall(directory)
    return json.loads((directory / "manifest.json").read_text()), sorted(p.name for p in (directory / "resources").iterdir())


def _run(directory, replacements=None):
    result = subprocess.run([sys.executable, "-W", "always", "-c", RUNNER, json.dumps(replacements or {})],
                            cwd=directory, capture_output=True, text=True, timeout=300)
    return result.stdout + result.stderr


def _with_path(path, optional):
    old = "data_path = ''" if optional else "Path('YOUR_ORIGINAL_DATA_FILE')"
    return {old: f"data_path = {str(path)!r}" if optional else f"Path({str(path)!r})"}


def _boxplot():
    from app.services.notebook_qc import render_qc
    return render_qc("protein", plot="boxplot", dpi=80)[1]


def _imputed_volcano():
    from app.schemas.stats import VolcanoRequest
    from app.services.notebook_volcano import render_volcano
    payload = VolcanoRequest(kind="protein", condition1="A", condition2="B", dataSource="imputed")
    return render_volcano(payload, width_cm=10, height_cm=8, dpi=80)[1]


def test_by_default_an_export_ships_no_data_matrix_and_needs_the_original_file(analysis, tmp_path):
    manifest, resources = _export(_boxplot(), tmp_path / "export")
    assert "input.parquet" not in resources and "table.parquet" in resources and "original.png" in resources
    assert manifest["provenance"]["requiredOriginalFiles"] == ["data_path"]
    assert set(manifest["input"]) >= {"sourceFingerprint", "inputFingerprint", "columns", "filename"}

    output = _run(tmp_path / "export")
    assert "Set data_path" in output and "ships no data matrix" in output

    output = _run(tmp_path / "export", _with_path(analysis, optional=False))
    assert "Fingerprint matches: the original file." in output
    assert "Every output is byte-identical to the one CoPYlot produced." in output


def test_a_shipped_matrix_runs_as_is_and_the_original_file_takes_over_when_named(analysis, tmp_path):
    manifest, resources = _export(_boxplot(), tmp_path / "export", ship_matrices=True)
    assert "input.parquet" in resources and manifest["provenance"]["requiredOriginalFiles"] == []
    notebook = next((tmp_path / "export").glob("*.ipynb")).read_text()
    assert "data_path = ''" in notebook

    shipped = _run(tmp_path / "export")
    assert "Fingerprint matches: the original file." not in shipped
    assert "Every output is byte-identical" in shipped

    named = _run(tmp_path / "export", _with_path(analysis, optional=True))
    assert "Fingerprint matches: the original file." in named and "Every output is byte-identical" in named


def test_another_file_is_reported_by_its_fingerprint_and_the_last_cell_decides(analysis, tmp_path):
    _export(_boxplot(), tmp_path / "export")
    other = pd.read_csv(analysis, sep="\t")
    other.loc[0, "A1"] = other.loc[0, "A1"] + 3.0
    other_path = tmp_path / "other.tsv"
    other.to_csv(other_path, sep="\t", index=False)
    output = _run(tmp_path / "export", _with_path(other_path, optional=False))
    assert "The original file is not the table CoPYlot used" in output
    assert "Some outputs differ" in output or "AssertionError" in output


@pytest.mark.parametrize("matrices, imputed, shipped, required", [
    (False, False, set(), ["data_path"]),
    (True, False, {"input_before_imputation.parquet"}, []),
    (False, True, {"input.parquet"}, []),
    (True, True, {"input.parquet"}, []),
])
def test_imputed_exports_ship_the_drawn_values_and_the_matrix_only_on_request(analysis, tmp_path, matrices, imputed, shipped, required):
    manifest, resources = _export(_imputed_volcano(), tmp_path / "export", ship_matrices=matrices, ship_imputed=imputed)
    assert "imputed_values.parquet" in resources
    assert {name for name in resources if name in {"input.parquet", "input_before_imputation.parquet"}} == shipped
    assert manifest["provenance"]["requiredOriginalFiles"] == required
    drawn = pd.read_parquet(tmp_path / "export" / "resources" / "imputed_values.parquet")
    assert list(drawn.columns) == ["row", "column", "value"] and len(drawn) > 0 and set(drawn["column"]) <= set(SAMPLES)

    if not required:
        as_shipped = _run(tmp_path / "export")
        assert "Every output is byte-identical" in as_shipped
        assert "Fill in the original-file path at the top to run the imputation itself" in as_shipped
    from_file = _run(tmp_path / "export", _with_path(analysis, optional=not required))
    assert "Every output is byte-identical" in from_file
    # The imputation call is shown and run; it draws exactly the shipped values here.
    assert "Imputing again here gives exactly the values CoPYlot drew." in from_file


def test_the_qq_plot_has_an_export(analysis, tmp_path):
    from app.services.notebook_distribution import render_qqnorm
    from app.services.plot_images import distribution_qqnorm_plot

    content, token = render_qqnorm("protein", random_seed=4)
    assert token and content == distribution_qqnorm_plot("protein", random_seed=4)
    manifest, resources = _export(token, tmp_path / "export")
    assert manifest["plot"] == "QQ Plot" and "input.parquet" not in resources
    assert "Every output is byte-identical" in _run(tmp_path / "export", _with_path(analysis, optional=False))


def test_viewers_get_no_export_until_the_one_switch_is_set(monkeypatch):
    from app.services import viewer_bootstrap, viewer_policy

    monkeypatch.setattr(viewer_bootstrap, "viewer_mode_enabled", lambda: True)
    assert viewer_policy.notebook_exports_blocked() is True
    assert not viewer_policy.viewer_request_allowed("GET", "/api/plots/notebooks/" + "a" * 32 + ".zip")
    monkeypatch.setattr(viewer_policy, "VIEWER_NOTEBOOK_EXPORTS", True)
    assert viewer_policy.notebook_exports_blocked() is False


def test_the_expected_file_name_stands_beside_the_path(analysis, tmp_path):
    _export(_boxplot(), tmp_path / "required")
    _export(_boxplot(), tmp_path / "optional", ship_matrices=True)
    first = lambda folder: json.loads(next((tmp_path / folder).glob("*.ipynb")).read_text())["cells"][1]["source"]  # noqa: E731
    assert 'data_path = Path(\'YOUR_ORIGINAL_DATA_FILE\')  # "protein.tsv" (required' in "".join(first("required"))
    assert 'data_path = \'\'  # "protein.tsv" (optional' in "".join(first("optional"))
    assert str(analysis.parent) not in "".join(first("required"))


@pytest.mark.parametrize("source", ["filtered", "imputed"])
def test_the_projection_has_an_export_with_its_table(analysis, tmp_path, source):
    from app.schemas.stats import SimulationRequest
    from app.services import stats_tools
    from app.services.notebook_projection import render_projection

    payload = SimulationRequest(kind="protein", condition1="A", condition2="B", dataSource=source,
                                sampleSizeOverride=6, varianceMultiplier=0.5)
    content, token = render_projection(payload, advanced={"titleElement": {"text": "What if"}})
    assert token and "What if" in content
    live = stats_tools._simulation_dataframe(payload)[2]
    manifest, resources = _export(token, tmp_path / "export")
    assert "input.parquet" not in resources and "table.parquet" in resources
    pd.testing.assert_frame_equal(pd.read_parquet(tmp_path / "export" / "resources" / "table.parquet"), live)
    output = _run(tmp_path / "export", _with_path(analysis, optional=False))
    assert "Every output is byte-identical" in output
    assert (tmp_path / "export" / "reproduced_table.csv").exists()


def test_the_kinase_library_volcano_has_an_export_of_the_uploaded_results(tmp_path):
    from app.services.notebook_kinase_volcano import render_kinase_volcano

    rng = np.random.default_rng(2)
    def results(prefix, count):
        return pd.DataFrame({"Kinase": [f"{prefix}{i}" for i in range(count)],
                             "Log2(Frequency Factor)": rng.normal(0, 1, count),
                             "Adjusted p-value": rng.uniform(0, 1, count) ** 3})
    content, token = render_kinase_volcano(st_results=results("STK", 30), tnc_results=results("TK", 10),
        st_name="st.csv", tnc_name="y.csv", p_value_threshold=0.1, header=True, condition1="A", condition2="B",
        condition1_control="", condition2_control="", annotate_hits=True, advanced=None)
    assert token and "copylot-plot" in content
    manifest, resources = _export(token, tmp_path / "export", ship_matrices=True, ship_imputed=True)
    assert {"st_results.parquet", "tnc_results.parquet", "original.html"} <= set(resources)
    assert manifest["provenance"]["uploadedFiles"] == {"st": "st.csv", "tnc": "y.csv"}
    assert "Every output is byte-identical" in _run(tmp_path / "export")


def test_a_peptide_export_follows_the_matrix_option(peptide_input, tmp_path):
    from app.services.notebook_peptide import render_coverage

    _, token = render_coverage(level="peptide", include_id=False, summary=False, text=True, text_size=12,
                               random_seed=3, quantity_column="Precursor.Quantity", dpi=80, advanced=None,
                               reference_line=False)
    manifest, resources = _export(token, tmp_path / "default")
    assert "input.parquet" not in resources and manifest["provenance"]["requiredOriginalFiles"] == ["data_path"]
    assert "ships no data matrix" in _run(tmp_path / "default")
    from_file = _run(tmp_path / "default", _with_path(peptide_input, optional=False))
    assert "Fingerprint matches: the original file." in from_file and "Every output is byte-identical" in from_file

    manifest, resources = _export(token, tmp_path / "shipped", ship_matrices=True)
    assert "input.parquet" in resources and manifest["provenance"]["requiredOriginalFiles"] == []
    assert "Every output is byte-identical" in _run(tmp_path / "shipped")
    assert "Every output is byte-identical" in _run(tmp_path / "shipped", _with_path(peptide_input, optional=True))


@pytest.fixture
def peptide_report(tmp_path, monkeypatch):
    from app.services import peptide_tools
    from app.services.dataset_store import save_peptide_path
    from app.services.peptide_metadata_store import save_peptide_metadata

    runs = ["ctrl_01.raw", "treated_01.raw"]
    stripped = ["AKPAK", "AKRAK", "PEPTIDEK", "MMK", "AAAK", "ABCK"]
    modified = ["AKPAK", "AKRAK", "(UniMod:1)PEPTIDEK", "M(UniMod:35)M(UniMod:35)K", "AAAK", "ABC(UniMod:4)K"]
    frame = pd.DataFrame({
        "File.Name": [runs[0]] * 6 + [runs[1]] * 6,
        "ProteinNames": ["TEST_HUMAN"] * 12,
        "Stripped.Sequence": stripped * 2,
        "Modified.Sequence": modified * 2,
        "Precursor.Quantity": [1., 3., 0., 5., 7., 10., 2., 0., 4., 2., 9., 10.],
        "RT": [11., 13., 15., 14., 17., 18., 11.5, 13.2, 15.1, 14.4, 17.3, 18.2],
        "Predicted.RT": [10., 12., 13., 15., 16., 20., 10.4, 12.1, 13.3, 15.2, 16.6, 19.],
    })
    path = tmp_path / "report.csv"
    frame.to_csv(path, index=False)
    save_peptide_path(str(path))
    save_peptide_metadata("metadata.csv", pd.DataFrame({"sample": runs, "condition": ["A", "B"]}))
    database = pd.DataFrame([("sp|P00001|TEST_HUMAN", "MMKAAAKABCKPEPTIDEKAKPAKAKRAK", "TEST_HUMAN", "P00001")],
                            columns=["full_header", "V1", "name", "accession"])
    monkeypatch.setattr(peptide_tools, "_load_fasta_dataframe", lambda species: database.copy())
    return path


def _peptide_token(output):
    from app.services import notebook_peptide
    from app.services.notebook_sequence import capture_sequence_coverage

    if output == "retention time":
        return notebook_peptide.render_rt(method="Scatter Plot", bins=20, dpi=80, add_line=True, header=True,
                                          sampling_seed=21, advanced=None)[1]
    if output == "interactive retention time":
        return notebook_peptide.render_interactive_rt(add_line=True, header=True, sampling_seed=32, advanced=None)[1]
    if output == "modifications":
        return notebook_peptide.render_modification(include_id=False, dpi=80, header=True, legend=True, advanced=None)[1]
    if output == "missed cleavages":
        return notebook_peptide.render_missed_cleavage(protease="trypsin", dpi=80, include_id=False, text=True,
                                                       header=True, text_size=11, advanced=None)[1]
    return capture_sequence_coverage(species="human", protein="TEST_HUMAN", chunk_size=30)["notebookToken"]


@pytest.mark.parametrize("output", ["retention time", "interactive retention time", "modifications",
                                    "missed cleavages", "sequence coverage"])
def test_every_peptide_export_runs_from_the_file_and_from_a_shipped_matrix(peptide_report, tmp_path, output):
    token = _peptide_token(output)
    manifest, resources = _export(token, tmp_path / "default")
    assert "input.parquet" not in resources and manifest["provenance"]["requiredOriginalFiles"] == ["data_path"]
    from_file = _run(tmp_path / "default", _with_path(peptide_report, optional=False))
    assert "Fingerprint matches: the original file." in from_file
    assert "the one CoPYlot produced" in from_file or "checked exactly above" in from_file, from_file[-600:]

    manifest, resources = _export(token, tmp_path / "shipped", ship_matrices=True)
    assert "input.parquet" in resources and manifest["provenance"]["requiredOriginalFiles"] == []
    for replacements in (None, _with_path(peptide_report, optional=True)):
        shipped = _run(tmp_path / "shipped", replacements)
        assert "the one CoPYlot produced" in shipped or "checked exactly above" in shipped, shipped[-600:]


@pytest.fixture
def phospho_package(tmp_path):
    rng = np.random.default_rng(11)
    paths = {}
    sites = pd.DataFrame(rng.normal(12, 1, (30, 6)), columns=SAMPLES)
    sites.iloc[:, :3] += np.repeat(rng.normal(0, 1, 10), 3)[:, None]
    sites["Protein_group"] = [f"P{i // 3}" for i in range(30)]
    sites["Gene_group"] = sites["Protein_group"]
    sites["PTM_Collapse_key"] = [f"P{i // 3}~G{i // 3}_S{i}_M1" for i in range(30)]
    sites["PTM_localization"] = [f"S{i}" for i in range(30)]
    sites["UPD_seq"] = ["AAAAAAASAAAAAAA"] * 30
    proteins = pd.DataFrame(rng.normal(20, 1, (8, 6)), columns=SAMPLES)
    proteins.insert(0, "ProteinNames", [f"PN{i}" for i in range(8)])
    proteins.insert(1, "ProteinAccessions", [f"P{i}" for i in range(8)])
    for kind, frame in (("phospho", sites), ("protein", proteins)):
        frame.iloc[2, 0] = frame.iloc[5, 4] = np.nan
        paths[kind] = tmp_path / f"{kind}.tsv"
        frame.to_csv(paths[kind], sep="\t", index=False)
        with paths[kind].open("rb") as stream:
            save_table_dataset(paths[kind].name, kind, read_dataframe(paths[kind].name, stream), source_path=str(paths[kind]))
        generate_annotation(AnnotationGenerateRequest(kind=kind, isLog2Transformed=True, conditions=CONDITIONS))
        data_tools.run_imputation(ImputationRunRequest(kind=kind))
    return paths


@pytest.mark.parametrize("source", ["filtered", "imputed"])
def test_the_protein_corrected_projection_has_an_export_that_needs_both_files(phospho_package, tmp_path, source):
    from app.schemas.stats import SimulationRequest
    from app.services import stats_tools
    from app.services.notebook_projection import render_projection

    payload = SimulationRequest(kind="phospho", condition1="A", condition2="B", dataSource=source,
                                proteinCorrection=True, sampleSizeOverride=6, varianceMultiplier=0.5)
    content, token = render_projection(payload, advanced=None)
    assert token
    live = stats_tools._simulation_dataframe(payload)[2]
    assert "protein_log2FC" in live.columns

    manifest, resources = _export(token, tmp_path / "default")
    assert not {"input.parquet", "protein_input.parquet"} & set(resources)
    assert manifest["provenance"]["requiredOriginalFiles"] == ["data_path", "protein_data_path"]
    pd.testing.assert_frame_equal(pd.read_parquet(tmp_path / "default" / "resources" / "table.parquet"), live)
    both = {"data_path = Path('YOUR_ORIGINAL_DATA_FILE')": f"data_path = Path({str(phospho_package['phospho'])!r})",
            "protein_data_path = Path('YOUR_ORIGINAL_PROTEIN_FILE')": f"protein_data_path = Path({str(phospho_package['protein'])!r})"}
    output = _run(tmp_path / "default", both)
    assert "Fingerprint matches: the original file." in output and "Fingerprint matches: the original protein file." in output
    assert "Every output is byte-identical" in output, output[-600:]

    manifest, resources = _export(token, tmp_path / "shipped", ship_matrices=True, ship_imputed=True)
    assert {"input.parquet", "protein_input.parquet"} <= set(resources)
    assert manifest["provenance"]["requiredOriginalFiles"] == []
    assert "Every output is byte-identical" in _run(tmp_path / "shipped")
