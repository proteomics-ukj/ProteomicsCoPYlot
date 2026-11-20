import io

import pytest
from fastapi import HTTPException, UploadFile

from app.api import routes_annotations, routes_datasets, routes_external


@pytest.mark.parametrize("route,label", [
    ("metadata", "Metadata file"),
    ("phosprot", "Phosphoprotein file"),
    ("dataset", "Dataset"),
    ("peptide", "Peptide dataset"),
    ("collapse", "Peptide dataset"),
])
@pytest.mark.parametrize("filename,status", [(None, 400), ("", 400), ("large.csv", 413)])
def test_upload_routes_reject_before_consuming_or_parsing(route, label, filename, status, monkeypatch):
    monkeypatch.setenv("COPYLOT_UPLOAD_MAX_BYTES", "4")
    stream = io.BytesIO(b"too large to accept")
    stream.seek(2)
    upload = UploadFile(filename=filename, file=stream)

    with pytest.raises(HTTPException) as raised:
        if route == "metadata":
            routes_annotations.upload_annotation_metadata(file=upload, kind="protein")
        elif route == "phosprot":
            routes_annotations.upload_phosprot_annotation(file=upload, is_log2_transformed="")
        elif route == "dataset":
            routes_datasets.upload_dataset(file=upload, kind="protein", package_name=None)
        elif route == "peptide":
            routes_datasets.upload_peptide_file(file=upload, package_name=None)
        else:
            routes_external.upload_peptide_collapse_route(
                routes_external.PeptideCollapseUpload(file=upload, cutoff=0.0)
            )

    assert raised.value.status_code == status
    assert raised.value.detail == (
        "Missing filename" if status == 400 else f"{label} exceeds the configured upload limit."
    )
    assert stream.tell() == 2
    assert not stream.closed
    assert isinstance(raised.value.__cause__, ValueError) if status == 413 else raised.value.__cause__ is None
