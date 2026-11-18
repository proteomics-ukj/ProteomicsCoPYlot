"""Shared request validation before an uploaded file is parsed or stored."""

from fastapi import HTTPException, UploadFile

from app.services.upload_limits import ensure_upload_size


def require_upload(file: UploadFile, *, label: str) -> str:
    if not file.filename:
        raise HTTPException(status_code=400, detail="Missing filename")
    try:
        ensure_upload_size(file, label=label)
    except ValueError as e:
        raise HTTPException(status_code=413, detail=str(e)) from e
    return file.filename
