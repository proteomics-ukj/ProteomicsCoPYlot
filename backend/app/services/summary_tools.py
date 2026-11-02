from __future__ import annotations


from app.schemas.annotation import AnnotationKind
from app.schemas.summary import SummaryMetadataTable, SummaryTablesResponse
from app.services.annotation_store import get_annotation
from app.services.runtime_cache import apply_cached_wrappers
from app.services.json_records import json_records
from app.services.dataset_store import get_current_dataset
from app.services.peptide_metadata_store import get_peptide_metadata


def _metadata_table(kind: AnnotationKind) -> SummaryMetadataTable:
    stored = get_annotation(kind)
    if stored is None or stored.metadata.empty:
        return SummaryMetadataTable(kind=kind, available=False)

    table = stored.metadata
    return SummaryMetadataTable(
        kind=kind,
        available=True,
        rows=len(table),
        columns=len(table.columns),
        columnNames=[str(column) for column in table.columns],
        table=json_records(table),
    )


def summary_tables() -> SummaryTablesResponse:
    frames = [stored.metadata for kind in ("protein", "phospho")
              if (stored := get_annotation(kind)) is not None]
    if get_current_dataset("protein") is None and (peptide := get_peptide_metadata()) is not None:
        frames.append(peptide.frame)
    samples = dict.fromkeys(str(sample) for frame in frames if "sample" in frame
                            for sample in frame["sample"].dropna())
    return SummaryTablesResponse(
        protein=_metadata_table("protein"),
        phospho=_metadata_table("phospho"),
        sampleNames=list(samples),
    )


apply_cached_wrappers(globals(), ["summary_tables"])
