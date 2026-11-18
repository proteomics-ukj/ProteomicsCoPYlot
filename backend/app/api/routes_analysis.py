from __future__ import annotations

import logging

from fastapi import APIRouter

from app.api.error_handling import api_errors
from app.schemas.analysis import AnalysisVolcanoRequest, AnalysisVolcanoResponse
from app.services.analysis_functions import analysis_volcano_data

LOGGER = logging.getLogger(__name__)


router = APIRouter(prefix="/api/analysis", tags=["analysis"])


@router.post("/volcano/data", response_model=AnalysisVolcanoResponse)
def analysis_volcano_data_route(payload: AnalysisVolcanoRequest) -> AnalysisVolcanoResponse:
    with api_errors("Failed to load analysis volcano data", logger=LOGGER):
        return analysis_volcano_data(payload)
