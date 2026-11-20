"""Lossless metadata-profile workbook exchange, separate from ordinary grouping uploads."""
from __future__ import annotations

import io

import pandas as pd
from openpyxl import Workbook, load_workbook

from app.schemas.annotation import AnnotationKind
from app.services.annotation_store import get_annotation, list_annotation_profiles
from app.services.metadata_upload_store import (
    get_active_profile, get_uploaded_metadata, list_auto_generated_profiles,
    list_profiles, MAX_MANUAL_METADATA_PROFILES,
)
from app.services.dataset_reader import _ensure_workbook_within_budget
from app.services.upload_limits import ensure_table_within_cell_budget

PROFILE_SHEET = "CoPYlot profiles v1"
PROFILE_COLUMNS = ["profile", "kind", "active"]
SAMPLE_COLUMNS = ["profile", "sample", "condition", "pair_id"]


def export_metadata_profiles(kind: AnnotationKind) -> bytes:
    source_kind = "phospho" if kind == "phosprot" else kind
    names = list(dict.fromkeys([*list_profiles(active_first=False), *list_annotation_profiles()]))
    records = []
    automatic = set(list_auto_generated_profiles())
    for name in names:
        if name in automatic:
            continue
        annotation = get_annotation(source_kind, metadata_profile_name=name)
        uploaded = get_uploaded_metadata(source_kind, profile_name=name)
        frame = annotation.metadata if annotation is not None else uploaded.frame if uploaded is not None else None
        if name == "Metadata" and frame is None and len(names) > 1:
            continue
        # Empty named profiles are included so their names survive the round trip.
        records.append((name, frame))
    if not any(frame is not None and not frame.empty for _, frame in records):
        raise ValueError("No metadata profiles are available to download for this level.")
    workbook = Workbook()
    profiles = workbook.active
    profiles.title = PROFILE_SHEET
    samples = workbook.create_sheet("Samples")
    profiles.append(PROFILE_COLUMNS)
    samples.append(SAMPLE_COLUMNS)
    active = get_active_profile()
    if active not in {name for name, _ in records}:
        active = records[0][0]
    for name, frame in records:
        profiles.append([name, source_kind, "yes" if name == active else "no"])
        if frame is not None:
            for row in frame.to_dict("records"):
                samples.append([name, *("" if pd.isna(row.get(key, "")) else str(row.get(key, ""))
                                        for key in SAMPLE_COLUMNS[1:])])
    for sheet in workbook:
        sheet.freeze_panes = "A2"
        sheet.auto_filter.ref = sheet.dimensions
        for row in sheet:
            for cell in row:
                # Labels beginning '=' are text, never spreadsheet formulas.
                cell.data_type = "s"
                cell.number_format = "@"
        for column in ("A", "B", "C", "D"):
            sheet.column_dimensions[column].width = 32
    output = io.BytesIO()
    workbook.save(output)
    return output.getvalue()


def read_metadata_profiles(stream, kind: str, data: pd.DataFrame):
    """Return None for ordinary XLSX files; validate a bundle fully before changing state."""
    _ensure_workbook_within_budget(stream)
    try:
        workbook = load_workbook(stream, read_only=True, data_only=False)
    except Exception as error:  # noqa: BLE001 - malformed uploads must produce a validation error
        stream.seek(0)
        raise ValueError("Could not read the Excel metadata file. Check that it is a valid .xlsx workbook.") from error
    try:
        if PROFILE_SHEET not in workbook.sheetnames:
            return None
        if "Samples" not in workbook.sheetnames:
            raise ValueError("Profile workbook is missing the Samples sheet.")
        def rows(sheet, columns):
            values = workbook[sheet].iter_rows(values_only=True)
            if list(next(values, ())) != columns:
                raise ValueError(f"Unexpected columns in the {sheet} sheet.")
            result = []
            for row_number, row in enumerate(values, start=2):
                ensure_table_within_cell_budget(row_number, len(columns), label="Metadata workbook")
                if not any(value is not None for value in row):
                    continue
                result.append(["" if value is None else str(value) for value in row])
            return result
        profiles = rows(PROFILE_SHEET, PROFILE_COLUMNS)
        samples = rows("Samples", SAMPLE_COLUMNS)
        if not profiles or len(profiles) > MAX_MANUAL_METADATA_PROFILES:
            raise ValueError("Profile workbook must contain between 1 and 5 original profiles.")
        names = [row[0] for row in profiles]
        if any(not name or name != name.strip() for name in names) or len(set(names)) != len(names):
            raise ValueError("Profile names must be nonempty, unique and have no surrounding whitespace.")
        if any(row[1] != kind for row in profiles):
            raise ValueError(f"This workbook does not contain {kind} metadata. Upload it for the matching level.")
        if any(row[2] not in {"yes", "no"} for row in profiles) or sum(row[2] == "yes" for row in profiles) != 1:
            raise ValueError("Exactly one profile must be marked active with 'yes'.")
        if any(row[0] not in names for row in samples):
            raise ValueError("Samples sheet refers to an unknown profile.")
        records = []
        for name in names:
            frame = pd.DataFrame([row[1:] for row in samples if row[0] == name], columns=SAMPLE_COLUMNS[1:])
            if frame['sample'].duplicated().any() or (frame[['sample', 'condition']] == '').any().any():
                raise ValueError(f"Profile '{name}' has duplicate samples or empty sample/condition values.")
            if not set(frame['sample']) <= set(map(str, data.columns)):
                raise ValueError(f"Profile '{name}' contains sample names missing from the loaded dataset.")
            records.append((name, frame))
        return records, next(row[0] for row in profiles if row[2] == "yes")
    finally:
        workbook.close()
        stream.seek(0)
