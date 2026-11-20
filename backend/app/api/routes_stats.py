from __future__ import annotations

import logging

from fastapi import APIRouter, Request, Response
from pydantic_core import to_json

from app.api.error_handling import api_errors
from app.schemas.annotation import AnnotationKind
from app.schemas.stats import (
    AnovaRequest,
    AnovaResultResponse,
    EnrichmentOrganism,
    EnrichmentRequest,
    EnrichmentResultResponse,
    ListEnrichmentRequest,
    ListEnrichmentResultResponse,
    PathwayOptionsResponse,
    SimulationRequest,
    SimulationResultResponse,
    StatisticalOptionsResponse,
    VolcanoControlRequest,
    VolcanoResultResponse,
    VolcanoRequest,
)
from app.services.stats_tools import (
    run_anova,
    pathway_options,
    run_enrichment_from_list,
    run_enrichment,
    run_simulation,
    run_volcano,
    run_volcano_control,
    statistical_options,
)

LOGGER = logging.getLogger(__name__)


router = APIRouter(prefix="/api/stats", tags=["stats"])


_COMPACT_STATS = "application/vnd.copylot.stats+json"


def _stats_response(result, request: Request):
    if request.headers.get("accept") != _COMPACT_STATS:
        return result

    table_fields = {"rows", "untestedRows", "clusterRows", "postHocRows"} & type(result).model_fields.keys()
    data = result.model_dump(exclude=table_fields)
    for key in table_fields:
        rows = getattr(result, key)
        data[key] = rows
        if not rows:
            continue
        columns = tuple(rows[0])
        if all(isinstance(row, dict) and tuple(row) == columns for row in rows):
            data[key] = {"columns": columns, "values": [list(row.values()) for row in rows]}
    return Response(to_json(data, inf_nan_mode="null"), media_type=_COMPACT_STATS, headers={"Vary": "Accept"})


@router.get("/options/{kind}", response_model=StatisticalOptionsResponse)
def stats_options_route(kind: AnnotationKind) -> StatisticalOptionsResponse:
    with api_errors("Failed to load statistical options", logger=LOGGER):
        return statistical_options(kind)


@router.post("/volcano/run", response_model=VolcanoResultResponse)
def volcano_route(payload: VolcanoRequest, request: Request) -> VolcanoResultResponse:
    with api_errors("Failed to run volcano analysis", logger=LOGGER):
        return _stats_response(run_volcano(payload), request)


@router.post("/volcano-control/run", response_model=VolcanoResultResponse)
def volcano_control_route(payload: VolcanoControlRequest, request: Request) -> VolcanoResultResponse:
    with api_errors("Failed to run control volcano analysis", logger=LOGGER):
        return _stats_response(run_volcano_control(payload), request)


@router.post("/anova/run", response_model=AnovaResultResponse)
def anova_route(payload: AnovaRequest, request: Request) -> AnovaResultResponse:
    with api_errors("Failed to run ANOVA analysis", logger=LOGGER):
        return _stats_response(run_anova(payload), request)


@router.post("/ora/run", response_model=EnrichmentResultResponse)
def ora_route(payload: EnrichmentRequest) -> EnrichmentResultResponse:
    with api_errors("Failed to run enrichment analysis", logger=LOGGER):
        return run_enrichment(payload)


@router.post("/ora/list/run", response_model=ListEnrichmentResultResponse)
def ora_list_route(payload: ListEnrichmentRequest) -> ListEnrichmentResultResponse:
    with api_errors("Failed to run list enrichment analysis", logger=LOGGER):
        from app.services.assistant_snapshots import store_snapshot
        result = run_enrichment_from_list(payload)
        from app.services.notebook_ora import capture_gene_list
        try:
            result.notebookToken = capture_gene_list(payload, result)
        except (ValueError, ModuleNotFoundError):
            LOGGER.exception("Gene list enrichment notebook unavailable")
        data = result.model_dump(exclude={"evidenceToken", "notebookToken", "terms", "genes", "backgroundGenes"})
        data.update(type="canvas_ora", querySize=len(result.genes), totalRows=len(result.terms),
                    terms=[term.model_dump() for term in result.terms[:10]],
                    settings=payload.model_dump(exclude={"genes", "backgroundGenes", "namedList", "inputSource"}))
        result.evidenceToken = store_snapshot(data)
        return result


@router.get("/pathways", response_model=PathwayOptionsResponse)
def pathway_options_route(organism: EnrichmentOrganism = "human") -> PathwayOptionsResponse:
    with api_errors("Failed to load pathway options", logger=LOGGER):
        return pathway_options(organism)


@router.post("/simulation/run", response_model=SimulationResultResponse)
def simulation_route(payload: SimulationRequest, request: Request) -> SimulationResultResponse:
    with api_errors("Failed to run simulation", logger=LOGGER):
        return _stats_response(run_simulation(payload), request)
