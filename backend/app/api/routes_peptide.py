from __future__ import annotations

import logging

import pandas as pd
from fastapi import APIRouter, HTTPException

from app.api.error_handling import api_errors
from app.schemas.peptide import (
    PeptideCoverageResponse,
    PeptideMetadataGenerateRequest,
    PeptideMetadataResponse,
    PeptideOverviewResponse,
    PeptideQuantityOptionsResponse,
    PeptideSpecies,
)
from app.services.peptide_metadata_store import (
    clear_peptide_metadata,
    get_peptide_metadata,
    save_peptide_metadata,
)
from app.services.peptide_tools import (
    _file_name_column,
    get_peptide_frame,
    peptide_overview,
)
from app.services.peptide_quantity import default_peptide_quantity_column
from app.services.json_records import json_records

LOGGER = logging.getLogger(__name__)


router = APIRouter(prefix="/api/peptide", tags=["peptide"])


def _metadata_response(stored) -> PeptideMetadataResponse:
    return PeptideMetadataResponse(
        filename=stored.filename,
        rows=len(stored.frame),
        columns=len(stored.frame.columns),
        createdAt=stored.created_at,
        preview=json_records(stored.frame, 30),
    )


@router.get("/overview", response_model=PeptideOverviewResponse)
def get_peptide_overview() -> PeptideOverviewResponse:
    with api_errors("Failed to load peptide overview", logger=LOGGER):
        return PeptideOverviewResponse(**peptide_overview())


@router.get("/quantity-options", response_model=PeptideQuantityOptionsResponse)
def get_peptide_quantity_options() -> PeptideQuantityOptionsResponse:
    with api_errors("Failed to load peptide abundance options", logger=LOGGER):
        frame = get_peptide_frame()
        default_column = default_peptide_quantity_column(frame)
        return PeptideQuantityOptionsResponse(
            columns=[default_column] if default_column else [],
            defaultColumn=default_column,
            requiresSelection=False,
        )


@router.post("/metadata/generate", response_model=PeptideMetadataResponse)
def generate_peptide_metadata(
    payload: PeptideMetadataGenerateRequest,
) -> PeptideMetadataResponse:
    with api_errors("Failed to load peptide samples", logger=LOGGER):
        frame = get_peptide_frame()
        file_col = _file_name_column(frame)

    sample_values = frame[file_col].dropna().astype(str).str.strip()
    samples = sample_values[sample_values != ""].drop_duplicates().tolist()
    available = set(samples)
    assigned: set[str] = set()
    rows: list[dict[str, str]] = []
    for condition in payload.conditions:
        condition_name = condition.name.strip()
        if not condition_name:
            continue
        for raw_sample in condition.columns:
            sample = str(raw_sample).strip()
            if sample not in available or sample in assigned:
                continue
            assigned.add(sample)
            rows.append({"sample": sample, "condition": condition_name})

    if not rows:
        raise HTTPException(
            status_code=400,
            detail="Assign at least one peptide file name to a condition.",
        )

    metadata = pd.DataFrame(rows, columns=["sample", "condition"])
    stored = save_peptide_metadata("Peptide Annotation (Data > Annotation)", metadata)
    return _metadata_response(stored)


@router.get("/metadata/current", response_model=PeptideMetadataResponse)
def get_current_peptide_metadata() -> PeptideMetadataResponse:
    stored = get_peptide_metadata()
    if stored is None:
        raise HTTPException(status_code=404, detail="No peptide metadata is currently stored.")
    return _metadata_response(stored)


@router.delete("/metadata/current")
def clear_current_peptide_metadata() -> dict[str, object]:
    clear_peptide_metadata(active_package_only=True)
    return {"cleared": True}


@router.get("/sequence-coverage", response_model=PeptideCoverageResponse)
def get_sequence_coverage(
    species: PeptideSpecies,
    protein: str,
    chunkSize: int = 100,
) -> PeptideCoverageResponse:
    with api_errors("Failed to calculate sequence coverage", logger=LOGGER):
        from app.services.notebook_sequence import capture_sequence_coverage
        return PeptideCoverageResponse(
            **capture_sequence_coverage(species=species, protein=protein, chunk_size=chunkSize)
        )
