"""Test setup for synthetic protein data already on the log2 scale."""

from app.schemas.annotation import AnnotationFilterConfig
from app.services.annotation_store import save_annotation


def save_log2_protein_annotation(frame, metadata):
    return save_annotation(
        kind="protein",
        source_data=frame,
        log2_data=frame,
        filtered_data=frame,
        metadata=metadata,
        is_log2_transformed=True,
        metadata_source="manual",
        filter_config=AnnotationFilterConfig(),
        auto_detected=False,
        warnings=[],
    )
