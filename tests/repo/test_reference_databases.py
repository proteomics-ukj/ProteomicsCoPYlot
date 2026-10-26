"""The bundled reference databases are what the logged upstream releases build.

`data/reference-sources/manifest.json` is the record: every source file and
every output in `data/db/` must match the checksum written there, the
documentation must quote the output checksums, and the build must be
deterministic (gzip without timestamps, sorted rows).
"""
from __future__ import annotations

import gzip
import csv
import importlib.util
import json
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "data" / "reference-sources" / "manifest.json"

spec = importlib.util.spec_from_file_location("build_reference_databases", ROOT / "tools" / "build_reference_databases.py")
build_tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_tool)


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def test_every_source_and_output_matches_the_manifest() -> None:
    build_tool.verify()  # raises SystemExit with the list of mismatches


def test_the_manifest_records_provenance_for_every_source(manifest) -> None:
    assert set(manifest["sources"]) == set(build_tool.SOURCE_FILES)
    for name, entry in manifest["sources"].items():
        for key in ("url", "release", "retrieved", "sha256", "licence", "licence_url", "citation", "bytes"):
            assert entry.get(key), f"{name} lacks {key}"
        assert entry["url"].startswith("https://")
        assert build_tool.source_path(name, entry["file"]).suffix in {".gz", ".fasta"}


def test_fasta_sources_are_stored_once_with_unchanged_provenance(manifest) -> None:
    fasta_sources = {name for name, entry in manifest["sources"].items() if entry["file"].endswith((".fasta.gz", ".fasta"))}
    assert fasta_sources == build_tool.DB_ONLY_SOURCES
    for name in fasta_sources:
        entry = manifest["sources"][name]
        path = build_tool.source_path(name, entry["file"])
        assert path == ROOT / "data" / "db" / entry["file"]
        assert not (ROOT / "data" / "reference-sources" / entry["file"]).exists()
        assert build_tool.sha256_of(path) == entry["sha256"] == manifest["outputs"][entry["file"]]["sha256"]
        assert manifest["build"]["sources_sha256"][name] == entry["sha256"]


def test_every_output_is_listed_and_documented(manifest) -> None:
    assert set(manifest["outputs"]) == set(build_tool.OUTPUTS)
    readme = (ROOT / "data" / "db" / "README.md").read_text(encoding="utf-8")
    for name, entry in manifest["outputs"].items():
        assert (ROOT / "data" / "db" / name).is_file()
        assert entry["sha256"] in readme, f"README does not quote the checksum of {name}"
    assert manifest["build"]["sources_sha256"] == {n: s["sha256"] for n, s in manifest["sources"].items()}


def test_cambridge_categories_distinguish_gpm_and_scientist_decisions() -> None:
    with (ROOT / "data/db/crap_ccp_metadata.csv").open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == len({r["cambridge_id"] for r in rows}) == 125
    assert sum(r["match_method"] == "exact_accession" for r in rows) == 111
    by_accession = {r["uniprot_accession"]: r for r in rows}
    assert by_accession["P02769"]["gpm_category_id"] == "1"  # bovine albumin
    assert by_accession["P13645"]["gpm_category_id"] == "2"  # human KRT10
    assert by_accession["P00004"]["gpm_category_id"] == "3"  # horse cytochrome c
    assert by_accession["P01375"]["gpm_category_id"] == "4"  # human TNF, UPS standard
    assert by_accession["P32503"]["gpm_category_id"] == "5"  # yeast virus capsid
    decisions = [r for r in rows if r["match_method"] == "scientific_decision"]
    assert len(decisions) == 14
    assert all(r["gpm_category"] for r in rows)
    assert sum(r["gpm_category_id"] == "2" for r in decisions) == 3
    assert sum(r["gpm_category_id"] == "1" for r in decisions) == 11
    for row in decisions:
        assert row["category_source"] == "scientific decision"
        assert row["category_source_version"] == ""
        assert row["gpm_entry_name"] == ""
    assert len([r for r in rows if r["uniprot_accession"] == "000000"]) == 2
    assert "decision_date" not in rows[0]
    assert "decision_rationale" not in rows[0]


def test_category_extraction_rejects_changed_source_layout() -> None:
    with pytest.raises(ValueError, match="structure/version changed"):
        build_tool.extract_gpm_categories(b"<html>temporarily unavailable</html>")


def test_cambridge_fetch_rejects_changed_bytes_before_replacing_snapshot(tmp_path, monkeypatch) -> None:
    target = tmp_path / "crap_ccp.fasta"
    target.write_bytes(b"existing snapshot")
    monkeypatch.setattr(build_tool, "SOURCES", tmp_path)
    monkeypatch.setattr(build_tool, "source_path", lambda *args: target)
    monkeypatch.setattr(build_tool, "_uniprot_release", lambda: ("unused", "unused"))
    monkeypatch.setattr(build_tool, "_go_release", lambda: "unused")
    monkeypatch.setattr(build_tool, "_get", lambda url: (b"changed upstream bytes", {}))
    with pytest.raises(ValueError, match="Cambridge v1.0 checksum mismatch"):
        build_tool.fetch(["cambridge_contaminants"])
    assert target.read_bytes() == b"existing snapshot"


def test_scientist_decisions_refuse_changed_record_identity(tmp_path, monkeypatch) -> None:
    fasta = (ROOT / "data/db/crap_ccp.fasta").read_text(encoding="utf-8")
    (tmp_path / "crap_ccp.fasta").write_text(
        fasta.replace("cRAP126|000000|ENDOP_GLUC", "cRAP126|000000|OTHER_REAGENT"),
        encoding="utf-8")
    monkeypatch.setattr(build_tool, "DB_SOURCES", tmp_path)
    monkeypatch.setattr(build_tool, "DB", tmp_path)
    with pytest.raises(ValueError, match="record identities changed"):
        build_tool.build_cambridge_metadata()


def test_nothing_from_the_sources_folder_is_shipped() -> None:
    """Executables and the image bundle data/db only; the sources stay in the repository."""
    for script in ("packaging/windows/build_exe_windows.ps1", "packaging/macos/build_all_macos_variants.sh", "Dockerfile"):
        text = (ROOT / script).read_text(encoding="utf-8")
        assert "reference-sources" not in text


def test_gzip_outputs_carry_no_timestamp(tmp_path) -> None:
    first, second = tmp_path / "a.gz", tmp_path / "b.gz"
    build_tool._write_gz_text(first, "x\ty\n")
    build_tool._write_gz_text(second, "x\ty\n")
    assert first.read_bytes() == second.read_bytes()
    with gzip.open(first, "rb") as handle:
        assert handle.read() == b"x\ty\n"


def test_verify_refuses_an_empty_manifest(monkeypatch) -> None:
    monkeypatch.setattr(build_tool, "_read_manifest", lambda: {"sources": {}, "outputs": {}})
    with pytest.raises(SystemExit, match="every expected source"):
        build_tool.verify()


def test_offline_rebuild_reproduces_frozen_outputs(tmp_path, monkeypatch, manifest) -> None:
    """Rebuild in isolation, without changing the shipped data or contacting upstream."""
    manifest_path = tmp_path / "manifest.json"
    shutil.copyfile(MANIFEST, manifest_path)
    monkeypatch.setattr(build_tool, "MANIFEST", manifest_path)
    monkeypatch.setattr(build_tool, "DB", tmp_path)
    monkeypatch.setattr(build_tool, "BUILD_LOG", tmp_path / "BUILD.log")
    monkeypatch.setattr(build_tool, "_get", lambda *_a, **_kw: pytest.fail("offline build tried to fetch"))
    build_tool.build()
    for name, entry in manifest["outputs"].items():
        assert build_tool.sha256_of(tmp_path / name) == entry["sha256"], name


@pytest.mark.parametrize(("value", "expected"), [
    ("ENST00000454207.1;ENSMUST000000001.2", ""),
    ("ENSG00000000001.4;ENST00000000001.1", "ENSG00000000001"),
    ("ENSMUSG00000000001.3;", "ENSMUSG00000000001"),
])
def test_ensembl_gene_fallback_never_returns_transcript_ids(value, expected) -> None:
    assert build_tool._ensembl_gene_ids(value) == expected


def test_shipped_gene_column_contains_only_gene_ids() -> None:
    import pandas as pd

    genes = pd.read_csv(ROOT / "data/db/BioNamesDB.txt.gz", sep="\t", usecols=["ensembl_gene_id"],
                        dtype=str, keep_default_na=False)["ensembl_gene_id"]
    assert not genes.str.contains(r"(?:ENST|ENSMUST)\d").any()


def test_go_build_filters_negation_species_and_nonprotein_objects(tmp_path, monkeypatch):
    build_tool._write_gz_text(tmp_path / "ontology.gz",
        "[Term]\nid: GO:1\nname: example\nnamespace: biological_process\nalt_id: GO:old\n")
    rows = []
    for symbol, relation, kind, taxon in [
        ("GENE", "involved_in", "protein", "taxon:9606"),
        ("NEGATED", "NOT|involved_in", "protein", "taxon:9606"),
        ("RNA", "involved_in", "lncRNA", "taxon:9606"),
        ("Mouse", "involved_in", "protein", "taxon:10090"),
        ("P1", "involved_in", "protein", "taxon:9606"),
        ("P2", "involved_in", "protein", "taxon:9606"),
    ]:
        rows.append("\t".join(["UniProtKB", symbol, symbol, relation, "GO:old", "PMID:1", "IEA",
                              "", "P", "", "", kind, taxon, "20261012", "GOA", "", ""]))
    # Named genes have a distinct accession; only P1/P2 need the symbol fallback.
    rows = [row.replace("UniProtKB\t" + name + "\t", "UniProtKB\tACC_" + name + "\t")
            if name not in {"P1", "P2"} else row
            for row, name in zip(rows, ["GENE", "NEGATED", "RNA", "Mouse", "P1", "P2"])]
    build_tool._write_gz_text(tmp_path / "annotations.gz", "\n".join(rows) + "\n")
    monkeypatch.setattr(build_tool, "SOURCES", tmp_path)
    monkeypatch.setattr(build_tool, "_read_manifest", lambda: {"sources": {
        "go_ontology": {"file": "ontology.gz"}, "go_annotations_human": {"file": "annotations.gz"}}})
    table = build_tool.build_go("human", {"P1": "RESOLVED"})
    assert table["term"].tolist() == ["GO:1"]
    assert table["genes"].tolist() == ["GENE;RESOLVED"]
    with pytest.raises(ValueError, match="Unsupported species"):
        build_tool.build_go("rat")


def test_go_annotations_propagate_over_is_a_and_part_of_only(tmp_path) -> None:
    obo = tmp_path / "go.obo.gz"
    build_tool._write_gz_text(obo, "\n".join([
        "format-version: 1.2", "",
        "[Term]", "id: GO:1", "name: root", "namespace: biological_process", "",
        "[Term]", "id: GO:2", "name: process", "namespace: biological_process", "is_a: GO:1 ! root", "",
        "[Term]", "id: GO:3", "name: part", "namespace: biological_process", "relationship: part_of GO:2 ! process", "",
        "[Term]", "id: GO:4", "name: regulator", "namespace: biological_process", "relationship: regulates GO:2 ! process", "is_a: GO:1 ! root", "",
        "[Term]", "id: GO:5", "name: old", "namespace: biological_process", "is_obsolete: true", "",
        "[Typedef]", "id: part_of", "",
    ]) + "\n")
    terms, parents, alt = build_tool._parse_obo(obo)
    assert set(terms) == {"GO:1", "GO:2", "GO:3", "GO:4"}
    cache: dict = {}
    assert build_tool._ancestors("GO:3", parents, cache) == {"GO:2", "GO:1"}
    assert build_tool._ancestors("GO:4", parents, cache) == {"GO:1"}, "regulates must not be followed"
