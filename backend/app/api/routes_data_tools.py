from __future__ import annotations

import logging


from fastapi import APIRouter, HTTPException
from fastapi.responses import Response

from app.schemas.annotation import AnnotationKind
from app.schemas.data_tools import (
    PhosphoLocalizationCutoffRequest,
    PhosphoLocalizationCutoffResponse,
    ConditionPaletteResponse,
    ConditionPaletteUpdateRequest,
    DistributionSummaryResponse,
    IdTranslationRequest,
    IdTranslationResponse,
    ImputationResultResponse,
    ImputationRunRequest,
)
from app.schemas.completeness import CompletenessTablesResponse
from app.services.completeness_tools import (
    completeness_tables,
)
from app.services.data_tools import (
    distribution_summary,
    export_id_translation,
    get_condition_palette_config,
    imputed_dataframe,
    update_condition_palette_config,
    run_id_translation,
    run_imputation,
)
from app.api.error_handling import api_errors
from app.api.download_headers import attachment_disposition

LOGGER = logging.getLogger(__name__)


router = APIRouter(prefix="/api/data-tools", tags=["data-tools"])


@router.post("/imputation/run", response_model=ImputationResultResponse)
def run_imputation_route(payload: ImputationRunRequest) -> ImputationResultResponse:
    with api_errors("Imputation failed", logger=LOGGER):
        return run_imputation(payload)


@router.post("/imputation/download")
def download_imputation_route(payload: ImputationRunRequest) -> Response:
    with api_errors("Imputation download failed", logger=LOGGER):
        frame = imputed_dataframe(payload)
        csv_bytes = frame.to_csv(index=False).encode("utf-8")
        filename = (
            f"imputed_{payload.kind}"
            f"_{payload.method}"
            f"_q{payload.qValue}"
            f"_std{payload.adjustStd}"
            f"_down{payload.downShift}"
            f"_width{payload.width}"
            f"_seed{payload.seed}"
            f"_samplewise{int(payload.sampleWise)}.csv"
        )
        return Response(
            content=csv_bytes,
            media_type="text/csv",
            headers={"Content-Disposition": attachment_disposition(filename)},
        )


@router.get("/distribution/{kind}", response_model=DistributionSummaryResponse)
def distribution_route(
    kind: AnnotationKind,
) -> DistributionSummaryResponse:
    with api_errors("Distribution analysis failed", logger=LOGGER):
        return distribution_summary(kind)


@router.get("/completeness/{kind}/tables", response_model=CompletenessTablesResponse)
def completeness_tables_route(
    kind: AnnotationKind,
    outlierThreshold: float = 50.0,
    includeId: bool = True,
    snapshot: str = "",
) -> CompletenessTablesResponse:
    if snapshot:
        from app.services.notebook_completeness import snapshot_tables
        return snapshot_tables(kind, snapshot, include_id=includeId, outlier_threshold=outlierThreshold)
    with api_errors("Completeness analysis failed", logger=LOGGER):
        return completeness_tables(
            kind=kind,
            outlier_threshold=outlierThreshold,
            include_id=includeId,
        )


@router.post("/id-translation/run", response_model=IdTranslationResponse)
def id_translation_route(payload: IdTranslationRequest) -> IdTranslationResponse:
    with api_errors("ID translation failed", logger=LOGGER):
        return run_id_translation(payload)


@router.post("/id-translation/download")
def download_id_translation_route(payload: IdTranslationRequest) -> Response:
    with api_errors("ID translation download failed", logger=LOGGER):
        filename, content, media_type = export_id_translation(payload)
        return Response(
            content=content,
            media_type=media_type,
            headers={"Content-Disposition": attachment_disposition(filename)},
        )


@router.get("/condition-colors/{kind}", response_model=ConditionPaletteResponse)
def get_condition_colors_route(kind: AnnotationKind) -> ConditionPaletteResponse:
    with api_errors("Failed to load condition colors", logger=LOGGER):
        return get_condition_palette_config(kind)


@router.post("/condition-colors/{kind}", response_model=ConditionPaletteResponse)
def set_condition_colors_route(
    kind: AnnotationKind,
    payload: ConditionPaletteUpdateRequest,
) -> ConditionPaletteResponse:
    with api_errors("Failed to save condition colors", logger=LOGGER):
        return update_condition_palette_config(kind, payload)


# Olsen et al. 2006: a site is "class I" at a localization probability of at
# least 0.75. The cutoff default stays 0; the page says
# how many sites below this line a lower cutoff admits.
CLASS_ONE_LOCALIZATION = 0.75


def _localization_cutoff_response() -> PhosphoLocalizationCutoffResponse:
    from app.services.annotation_store import get_annotation
    from app.services.dataset_store import get_current_dataset
    from app.services.functions import LOCALIZATION_COLUMN, filter_localization
    from app.services.phospho_localization import get_localization_cutoff

    annotation = get_annotation("phospho")
    frame = getattr(annotation, "source_data", None) if annotation is not None else None
    if frame is None or getattr(frame, "empty", True):
        dataset = get_current_dataset("phospho")
        frame = getattr(dataset, "frame", None)
    available = frame is not None and LOCALIZATION_COLUMN in frame.columns
    cutoff = get_localization_cutoff()
    below, sites = 0, 0
    if available and frame is not None:
        import pandas as pd

        retained = filter_localization(frame, cutoff)
        probability = pd.to_numeric(retained[LOCALIZATION_COLUMN], errors="coerce")
        probability = probability.where(probability.between(0.0, 1.0))
        sites = len(retained)
        below = int((probability < CLASS_ONE_LOCALIZATION).sum())
    return PhosphoLocalizationCutoffResponse(
        cutoff=cutoff, available=bool(available), belowClassOne=below, sites=sites,
    )


@router.get("/phospho-localization-cutoff", response_model=PhosphoLocalizationCutoffResponse)
def get_phospho_localization_cutoff_route() -> PhosphoLocalizationCutoffResponse:
    return _localization_cutoff_response()


@router.post("/phospho-localization-cutoff", response_model=PhosphoLocalizationCutoffResponse)
def set_phospho_localization_cutoff_route(
    payload: PhosphoLocalizationCutoffRequest,
) -> PhosphoLocalizationCutoffResponse:
    from app.services.phospho_localization import set_localization_cutoff

    try:
        set_localization_cutoff(payload.cutoff)
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e)) from e
    return _localization_cutoff_response()
