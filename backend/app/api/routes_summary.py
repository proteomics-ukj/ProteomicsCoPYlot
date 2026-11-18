from __future__ import annotations

import logging

from fastapi import APIRouter, HTTPException

from app.api.error_handling import api_errors
from app.schemas.summary import (
    SummaryReportRequest,
    SummaryReportResponse,
    SummaryTablesResponse,
)
from app.services.plot_export import RENDER_LOCK
from app.services.report_function import report_function
from app.services.notebook_report import report_with_notebook
from app.services.summary_tools import summary_tables

LOGGER = logging.getLogger(__name__)


router = APIRouter(prefix="/api/summary", tags=["summary"])


@router.get("/tables", response_model=SummaryTablesResponse)
def summary_tables_route() -> SummaryTablesResponse:
    try:
        return summary_tables()
    except Exception as e:
        LOGGER.exception("Failed to load summary tables")
        raise HTTPException(status_code=500, detail="Failed to load summary tables.") from e


@router.post("/report", response_model=SummaryReportResponse)
def summary_report_route(payload: SummaryReportRequest) -> SummaryReportResponse:
    with api_errors("Failed to build report", logger=LOGGER):
        # A report draws many figures; they are drawn one at a time like any
        # other plot (RENDER_LOCK in plot_export.py).
        with RENDER_LOCK:
            return report_with_notebook(report_function, payload)
