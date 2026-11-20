"""A save file holds what cannot be recomputed, and nothing else.

Since 1.0.0b2 a save file carries only the values imputation drew, a translated
table is rebuilt from its file plus the recorded translations, and a collapsed
phosphoprotein table is collapsed again from the restored phosphosites. Each is
checked against a fingerprint of the saved table; a table that would not come
back identical is stored whole, as before. The limits on save files exist only
when an operator sets them.
"""

from __future__ import annotations

import gzip
import io
import json

import numpy as np
import pandas as pd
import pytest
from fastapi import UploadFile

from app.api.routes_annotations import generate_annotation
from app.schemas.annotation import AnnotationFilterConfig, AnnotationGenerateRequest
from app.schemas.data_tools import IdTranslationRequest, ImputationRunRequest
from app.services import annotation_store, data_tools, id_translation, phosprot_tools, restore_state
from app.services.annotation_processor import UPLOADED_METADATA_NOTICE
from app.services.dataset_reader import read_dataframe
from app.services.dataset_store import get_current_dataset, reset_all_datasets, save_table_dataset
from app.services.metadata_upload_store import reset_all_uploaded_metadata

SAMPLES = ["A1", "A2", "A3", "B1", "B2", "B3"]
CONDITIONS = [{"name": group, "columns": [f"{group}{i}" for i in range(1, 4)]} for group in "AB"]


@pytest.fixture(autouse=True)
def isolated_stores():
    def reset():
        annotation_store.reset_all_annotations()
        reset_all_datasets()
        reset_all_uploaded_metadata()
    reset()
    yield
    reset()


def _export(portable: bool = False) -> dict:
    return json.loads(gzip.decompress(restore_state.export_state_blob(portable=portable)).decode("utf-8"))


def _clear() -> None:
    annotation_store.reset_all_annotations()
    reset_all_datasets()


def _load_protein(tmp_path, *, missing: int = 8):
    values = np.random.default_rng(7).uniform(8, 16, (12, 6))
    values.flat[np.random.default_rng(8).choice(values.size, missing, replace=False)] = np.nan
    raw = pd.DataFrame(values, columns=SAMPLES)
    raw["ProteinNames"] = [f"P{i}" for i in range(12)]
    path = tmp_path / "protein.tsv"
    raw.to_csv(path, sep="\t", index=False)
    with path.open("rb") as stream:
        frame = read_dataframe(path.name, stream)
    save_table_dataset(path.name, "protein", frame, source_path=str(path))
    return path


def _annotate(kind: str = "protein") -> None:
    generate_annotation(AnnotationGenerateRequest(kind=kind, isLog2Transformed=True, conditions=CONDITIONS))


def _replacement(path) -> restore_state.ReplacementInput:
    return restore_state.ReplacementInput(file=UploadFile(file=io.BytesIO(path.read_bytes()), filename=path.name))


# --- format -----------------------------------------------------------------


def test_a_save_file_names_its_format_by_the_release_version(tmp_path):
    _load_protein(tmp_path)
    state = _export()
    assert state["format"] == "copylot_state"
    assert state["formatVersion"] == restore_state._STATE_FORMAT_VERSION == "1.0.0b2"


@pytest.mark.parametrize("version, readable", [
    ("1.0.0b1", True), ("1.0.0b2", True), ("1.0.0-beta.2", True), ("v1.0.0-beta.1", True),
    ("1.0.0b3", False), ("1.0.0", False), ("1.1.0", False), ("2.0.0b1", False), ("", False), ("next", False),
])
def test_a_newer_format_is_refused_and_an_older_one_is_read(version, readable):
    blob = gzip.compress(json.dumps({"format": "copylot_state", "formatVersion": version}).encode())
    if readable:
        assert restore_state.decode_restore_blob(blob)["formatVersion"] == version
        return
    with pytest.raises(restore_state.RestoreStateError) as failure:
        restore_state.decode_restore_blob(blob)
    if version in {"1.0.0b3", "1.0.0", "1.1.0", "2.0.0b1"}:
        assert f"needs app version {version} or newer" in str(failure.value)


def test_releases_are_ordered_as_in_the_release_notes():
    ordered = ["1.0.0b1", "1.0.0b2", "1.0.0b10", "1.0.0", "1.0.1", "1.1.0", "2.0.0b1", "2.0.0"]
    keys = [restore_state._version_key(version) for version in ordered]
    assert keys == sorted(keys) and len(set(keys)) == len(keys)


def test_a_save_file_from_before_the_format_version_is_still_read(tmp_path):
    _load_protein(tmp_path)
    _annotate()
    data_tools.run_imputation(ImputationRunRequest(kind="protein"))
    saved = annotation_store.get_annotation("protein").imputed_data.copy()
    state = _export()
    # As 1.0.0b1 wrote it: the old format name, no version, the whole matrix.
    state["format"] = "copylot_state_v1"
    del state["formatVersion"]
    entry = state["backend"]["packages"][0]["annotations"][0]
    del entry["imputedSparse"]
    entry["imputedFrame"] = restore_state._frame_to_payload(saved)
    state = restore_state.decode_restore_blob(gzip.compress(json.dumps(state).encode()))
    _clear()
    result = restore_state.apply_restore_state(state, {})
    assert not any("Discarded" in warning for warning in result.warnings), result.warnings
    restored = annotation_store.get_annotation("protein")
    assert restored.imputation.applied is True
    np.testing.assert_allclose(restored.imputed_data[SAMPLES].to_numpy(float), saved[SAMPLES].to_numpy(float))


# --- limits -----------------------------------------------------------------


_LIMITS = [
    "COPYLOT_RESTORE_UPLOAD_MAX_BYTES", "COPYLOT_RESTORE_EXPANDED_MAX_BYTES",
    "COPYLOT_RESTORE_REPLACEMENT_MAX_BYTES", "COPYLOT_RESTORE_REPLACEMENTS_TOTAL_MAX_BYTES",
    "COPYLOT_RESTORE_TABLE_MAX_ROWS", "COPYLOT_RESTORE_TABLE_MAX_COLUMNS", "COPYLOT_RESTORE_TABLE_MAX_CELLS",
    "COPYLOT_RESTORE_REPLACEMENT_MAX_FILES", "COPYLOT_RESTORE_MANIFEST_MAX_CHARS",
]


def test_nothing_is_limited_until_an_operator_sets_a_limit(monkeypatch):
    for setting in _LIMITS:
        monkeypatch.delenv(setting, raising=False)
    # Far beyond every limit this used to have built in.
    restore_state._validate_table_dimensions(10**9, 10**6)
    files = [UploadFile(file=io.BytesIO(b"x"), filename=f"{index}.tsv") for index in range(300)]
    manifest = json.dumps([{"packageName": f"P{index}", "kind": "protein"} for index in range(300)])
    assert len(restore_state.parse_replacements_manifest(manifest + " " * 2_000_000, files)) == 300
    payload = {"format": "copylot_state", "formatVersion": "1.0.0b2", "uiState": {"text": "x" * 1000}}
    blob = restore_state._serialize_state_payload(payload)
    assert restore_state.decode_restore_blob(blob) == payload

    monkeypatch.setenv("COPYLOT_RESTORE_TABLE_MAX_ROWS", "5")
    with pytest.raises(restore_state.RestoreStateError, match="COPYLOT_RESTORE_TABLE_MAX_ROWS"):
        restore_state._validate_table_dimensions(6, 1)
    monkeypatch.setenv("COPYLOT_RESTORE_REPLACEMENT_MAX_FILES", "299")
    with pytest.raises(restore_state.RestoreStateError, match="Too many replacement files"):
        restore_state.parse_replacements_manifest(manifest, files)
    monkeypatch.setenv("COPYLOT_RESTORE_UPLOAD_MAX_BYTES", str(len(blob) - 1))
    with pytest.raises(restore_state.RestoreStateError, match="upload limit"):
        restore_state.decode_restore_blob(blob)


def test_a_replacement_file_is_not_bound_by_the_upload_limit(tmp_path, monkeypatch):
    path = _load_protein(tmp_path)
    state = _export(portable=True)
    _clear()
    monkeypatch.setenv("COPYLOT_UPLOAD_MAX_BYTES", "10")
    restore_state.apply_restore_state(state, {("Data", "protein"): _replacement(path)}, allow_local_paths=False)
    assert get_current_dataset("protein") is not None


# --- the notice of an uploaded metadata table --------------------------------


@pytest.mark.parametrize("source", ["manual", "auto", "uploaded"])
def test_a_restored_annotation_only_names_an_uploaded_table_when_one_was_used(tmp_path, source):
    _load_protein(tmp_path)
    _annotate()
    state = _export()
    entry = state["backend"]["packages"][0]["annotations"][0]
    entry["metadataSource"] = source
    del entry["warnings"]  # a save file from before the notices were saved
    _clear()
    restore_state.apply_restore_state(state, {})
    restored = annotation_store.get_annotation("protein")
    assert restored.metadata_source == source
    assert (UPLOADED_METADATA_NOTICE in restored.warnings) is (source == "uploaded")


def test_the_notices_of_an_annotation_come_back_as_they_were(tmp_path):
    _load_phospho(tmp_path, imputed=False)
    phosprot_tools.aggregate_from_phospho("mean")
    saved = {kind: list(annotation_store.get_annotation(kind).warnings) for kind in ("phospho", "phosprot")}
    assert any("aggregated from" in warning for warning in saved["phosprot"])
    state = _export()
    _clear()
    restore_state.apply_restore_state(state, {})
    assert {kind: annotation_store.get_annotation(kind).warnings for kind in saved} == saved


# --- imputed values -----------------------------------------------------------


def test_sparse_restore_rejects_changed_cells_even_when_column_sums_agree():
    basis = pd.DataFrame({"ProteinNames": ["P1", "P2", "P3"], "A": [10., 20., np.nan]})
    imputed = basis.fillna({"A": 15.})
    payload = restore_state._sparse_imputation_payload(imputed, basis, ["A"])
    changed = basis.copy()
    changed.loc[0, "A"], changed.loc[1, "A"] = 20., 10.
    restored, problem = restore_state._imputed_from_sparse(payload, changed)
    assert restored is None and problem is not None


def test_only_the_drawn_values_are_saved_and_the_matrix_comes_back_identical(tmp_path):
    _load_protein(tmp_path)
    _annotate()
    result = data_tools.run_imputation(ImputationRunRequest(kind="protein"))
    saved = annotation_store.get_annotation("protein").imputed_data.copy()
    entry = _export()["backend"]["packages"][0]["annotations"][0]
    assert entry["imputedFrame"] is None
    sparse = entry["imputedSparse"]
    drawn = sum(len(cells["values"]) for cells in sparse["cells"].values())
    assert drawn == result.missingBefore - result.missingAfter > 0
    assert set(sparse["cells"]) <= set(SAMPLES)

    state = _export()
    _clear()
    restored_result = restore_state.apply_restore_state(state, {})
    assert not any("Discarded" in warning for warning in restored_result.warnings), restored_result.warnings
    restored = annotation_store.get_annotation("protein")
    assert restored.imputation.applied is True
    pd.testing.assert_frame_equal(restored.imputed_data, saved)
    assert restore_state._table_fingerprint(restored.imputed_data) == sparse["sha256"]


def test_a_matrix_that_is_not_its_basis_plus_drawn_values_is_saved_whole(tmp_path):
    _load_protein(tmp_path)
    _annotate()
    data_tools.run_imputation(ImputationRunRequest(kind="protein"))
    annotation = annotation_store.get_annotation("protein")
    annotation.imputed_data = annotation.imputed_data.assign(Extra="x")
    entry = _export()["backend"]["packages"][0]["annotations"][0]
    assert "imputedSparse" not in entry
    assert "Extra" in entry["imputedFrame"]["columns"]


@pytest.mark.parametrize("sparse", [False, True])
@pytest.mark.parametrize("shift", [1e-12, 0.01])
def test_drawn_values_are_kept_only_on_the_values_they_were_drawn_for(tmp_path, shift, sparse):
    path = _load_protein(tmp_path)
    _annotate()
    data_tools.run_imputation(ImputationRunRequest(kind="protein"))
    state = _export()
    if not sparse:
        entry = state["backend"]["packages"][0]["annotations"][0]
        entry["imputedFrame"] = restore_state._frame_to_payload(annotation_store.get_annotation("protein").imputed_data)
        del entry["imputedSparse"]
    kept = not sparse and shift < 1e-9
    # The file is the same one (same fingerprint in the save file), read back
    # with observed values that moved. A full save supports cell-wise tolerance;
    # a sparse save cannot verify that tolerance from column totals alone.
    for entry in state["backend"]["packages"][0]["datasets"].values():
        if entry:
            for key in ("sizeBytes", "sha256", "tableSha256"):
                entry.pop(key, None)
    raw = pd.read_csv(path, sep="\t")
    raw[SAMPLES] = raw[SAMPLES] + shift
    raw.to_csv(path, sep="\t", index=False, float_format="%.17g")
    _clear()
    result = restore_state.apply_restore_state(state, {})
    restored = annotation_store.get_annotation("protein")
    assert restored.imputation.applied is kept
    assert restored.imputed_data.empty is (not kept)
    assert any("observed values" in warning or "matrix fingerprint" in warning for warning in result.warnings) is (not kept)


# --- a translated table -------------------------------------------------------


def _translate(monkeypatch) -> None:
    database = pd.DataFrame({
        "uniprot_symbol": [f"P{i}" for i in range(12)], "symbol": [f"G{i}" for i in range(12)],
        "ensembl_gene_id": [f"E{i}" for i in range(12)], "uniprot_gene": [f"UPG{i}" for i in range(12)],
    })
    monkeypatch.setattr(id_translation, "load_db", lambda: database)
    id_translation.run_id_translation(IdTranslationRequest(
        kind="protein", column="ProteinNames", inputDb="UniProt Accession", outputDb="Ensembl", autoDetectInput=False,
    ))


@pytest.mark.parametrize("portable", [False, True])
def test_a_translated_table_is_rebuilt_from_its_file(tmp_path, monkeypatch, portable):
    path = _load_protein(tmp_path)
    _translate(monkeypatch)
    translated = get_current_dataset("protein")
    assert translated.source_path == str(path) and "EnsemblNames" in translated.frame.columns
    saved = translated.frame.copy()
    history = [step.model_dump() for step in translated.id_translations]

    state = _export(portable=portable)
    entry = state["backend"]["packages"][0]["datasets"]["protein"]
    assert "frame" not in entry and "derived" not in entry
    assert entry["path"] == (None if portable else str(path))
    _clear()
    monkeypatch.setattr(id_translation, "load_db", lambda: pytest.fail("A restore must not query today's database"))
    if portable:
        missing = restore_state.inspect_missing_files(state, {}, allow_local_paths=False)
        assert [(item.package_name, item.kind) for item in missing] == [("Data", "protein")]
    replacements = {("Data", "protein"): _replacement(path)} if portable else {}
    result = restore_state.apply_restore_state(state, replacements, allow_local_paths=not portable)
    assert not any("differs" in warning or "translation" in warning for warning in result.warnings), result.warnings
    restored = get_current_dataset("protein")
    pd.testing.assert_frame_equal(restored.frame, saved)
    assert [step.model_dump() for step in restored.id_translations] == history


def test_a_translated_table_its_file_does_not_reproduce_is_saved_whole(tmp_path, monkeypatch):
    path = _load_protein(tmp_path)
    _translate(monkeypatch)
    saved = get_current_dataset("protein").frame.copy()
    path.write_text(path.read_text().replace("P3", "P300"))
    entry = _export(portable=True)["backend"]["packages"][0]["datasets"]["protein"]
    assert entry["derived"] is True and entry["frame"]["columns"] == list(saved.columns)
    state = _export(portable=True)
    _clear()
    assert restore_state.inspect_missing_files(state, {}, allow_local_paths=False) == []
    restore_state.apply_restore_state(state, {}, allow_local_paths=False)
    assert restore_state._table_fingerprint(get_current_dataset("protein").frame) == restore_state._table_fingerprint(saved)


def test_a_file_the_translation_no_longer_fits_is_restored_without_it(tmp_path, monkeypatch):
    path = _load_protein(tmp_path)
    _translate(monkeypatch)
    state = _export()
    path.write_text(path.read_text().replace("P3", "P300"))
    _clear()
    result = restore_state.apply_restore_state(state, {})
    assert any("identifier translation" in warning and "without it" in warning for warning in result.warnings)
    restored = get_current_dataset("protein")
    assert restored.id_translations == [] and "EnsemblNames" not in restored.frame.columns


# --- a collapsed phosphoprotein table ------------------------------------------


def _load_phospho(tmp_path, *, imputed: bool):
    values = np.random.default_rng(11).uniform(8, 16, (10, 6))
    values[0, 0] = values[3, 4] = np.nan
    sites = pd.DataFrame(values, columns=SAMPLES)
    sites.insert(0, "PTM_Collapse_key", [f"K{i}" for i in range(10)])
    sites.insert(0, "Gene_group", [f"G{i // 2}" for i in range(10)])
    sites.insert(0, "Protein_group", [f"P{i // 2}" for i in range(10)])
    sites.insert(0, "PTM_localization", ["S1"] * 10)
    sites.insert(0, "UPD_seq", [f"SEQ{i}" for i in range(10)])
    path = tmp_path / "sites.tsv"
    sites.to_csv(path, sep="\t", index=False)
    with path.open("rb") as stream:
        frame = read_dataframe(path.name, stream)
    save_table_dataset(path.name, "phospho", frame, source_path=str(path))
    annotation_store.save_annotation(
        kind="phospho", source_data=frame, metadata=pd.DataFrame({"sample": SAMPLES, "condition": list("AAABBB")}),
        log2_data=frame, filtered_data=frame, is_log2_transformed=True, metadata_source="manual",
        filter_config=AnnotationFilterConfig(minPresent=0), auto_detected=False, warnings=[],
    )
    if imputed:
        data_tools.run_imputation(ImputationRunRequest(kind="phospho"))
    return path


@pytest.mark.parametrize("mode", ["robust_site_median", "mean", "sum_ignore_na", "sum_propagate_na"])
@pytest.mark.parametrize("source", ["non_imputed", "imputed"])
@pytest.mark.parametrize("portable", [False, True])
def test_a_collapsed_table_is_collapsed_again_on_restore(tmp_path, mode, source, portable):
    path = _load_phospho(tmp_path, imputed=source == "imputed")
    phosprot_tools.aggregate_from_phospho(mode, source=source)
    saved = get_current_dataset("phosprot")
    frame, lineage = saved.frame.copy(), saved.phosprot_lineage.model_dump()
    filtered = annotation_store.get_annotation("phosprot").filtered_data.copy()

    state = _export(portable=portable)
    entry = state["backend"]["packages"][0]["datasets"]["phosprot"]
    assert "frame" not in entry and entry["path"] is None
    assert entry["derived"] is True and entry["recompute"]["type"] == "phosphoCollapse"
    _clear()
    missing = restore_state.inspect_missing_files(state, {}, allow_local_paths=not portable)
    assert [item.kind for item in missing] == (["phospho"] if portable else [])
    replacements = {("Data", "phospho"): _replacement(path)} if portable else {}
    result = restore_state.apply_restore_state(state, replacements, allow_local_paths=not portable)
    assert not any("rebuilt" in warning or "Discarded" in warning for warning in result.warnings), result.warnings
    assert "Data:phosprot" in result.restored_datasets
    restored = get_current_dataset("phosprot")
    pd.testing.assert_frame_equal(restored.frame, frame)
    assert restored.phosprot_lineage.model_dump() == lineage and restored.phosprot_aggregation == mode
    annotation = annotation_store.get_annotation("phosprot")
    assert annotation is not None and annotation.phosprot_aggregation == mode
    pd.testing.assert_frame_equal(annotation.filtered_data, filtered)


def test_a_collapsed_table_the_phosphosites_no_longer_give_is_saved_whole(tmp_path):
    _load_phospho(tmp_path, imputed=False)
    phosprot_tools.aggregate_from_phospho("mean")
    saved = get_current_dataset("phosprot").frame.copy()
    # The phosphosites were filtered differently after the collapse.
    phospho = annotation_store.get_annotation("phospho")
    phospho.filtered_data = phospho.filtered_data.iloc[:4].copy()
    entry = _export(portable=True)["backend"]["packages"][0]["datasets"]["phosprot"]
    assert "recompute" not in entry and entry["derived"] is True
    assert entry["frame"]["columns"] == list(saved.columns)


def test_a_collapse_that_cannot_be_repeated_is_reported(tmp_path):
    _load_phospho(tmp_path, imputed=True)
    phosprot_tools.aggregate_from_phospho("mean", source="imputed")
    state = _export()
    # The drawn phosphosite values are gone, so the imputed collapse cannot be repeated.
    for annotation in state["backend"]["packages"][0]["annotations"]:
        if annotation["kind"] == "phospho":
            annotation["imputedSparse"] = None
            annotation["imputation"] = None
    _clear()
    result = restore_state.apply_restore_state(state, {})
    assert any("could not be rebuilt" in warning and "Collapse the phosphosites again" in warning
               for warning in result.warnings), result.warnings
    assert get_current_dataset("phosprot") is None
    assert get_current_dataset("phospho") is not None


# --- saving again after a restore ---------------------------------------------


@pytest.mark.parametrize("portable", [False, True])
def test_saving_again_after_a_restore_gives_the_same_small_file(tmp_path, portable):
    path = _load_phospho(tmp_path, imputed=True)
    phosprot_tools.aggregate_from_phospho("mean", source="imputed")
    state = _export(portable=portable)
    _clear()
    replacements = {("Data", "phospho"): _replacement(path)} if portable else {}
    restore_state.apply_restore_state(state, replacements, allow_local_paths=not portable)
    restored = get_current_dataset("phospho")
    if portable:
        # The uploaded replacement is the dataset's file from here on, as an upload is.
        assert restored.source_path and restored.source_path != str(path)
    again = _export(portable=portable)
    for saved in (state, again):
        datasets = saved["backend"]["packages"][0]["datasets"]
        assert "frame" not in datasets["phospho"] and "frame" not in datasets["phosprot"]
        assert all(entry["imputedFrame"] is None for entry in saved["backend"]["packages"][0]["annotations"])
    assert again["backend"] == state["backend"] or not portable
    if portable:
        assert again["backend"]["packages"][0]["datasets"] == state["backend"]["packages"][0]["datasets"]
        assert again["backend"]["packages"][0]["annotations"] == state["backend"]["packages"][0]["annotations"]
