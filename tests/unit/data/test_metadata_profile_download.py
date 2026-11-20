"""Download/upload round trips preserve original profile names and per-profile pairing."""
import io

import pandas as pd
import pytest
from fastapi.testclient import TestClient
from openpyxl import load_workbook

from app.main import app
from app.services.annotation_store import reset_all_annotations, save_annotation
from app.schemas.annotation import AnnotationFilterConfig
from app.services.dataset_store import reset_all_datasets, save_table_dataset
from app.services.metadata_profile_io import export_metadata_profiles, read_metadata_profiles, PROFILE_SHEET
from app.services.metadata_upload_store import (
    get_active_profile, get_uploaded_metadata, list_auto_generated_profiles,
    reset_all_uploaded_metadata, save_uploaded_metadata_profiles, set_active_profile,
)


@pytest.fixture(autouse=True)
def clean_state():
    reset_all_datasets()
    reset_all_uploaded_metadata()
    reset_all_annotations()
    yield
    reset_all_datasets()
    reset_all_uploaded_metadata()
    reset_all_annotations()


def setup_profiles(kind='protein', count=2):
    data = pd.DataFrame({'ProteinNames': ['P1', 'P2'], '001': [1., 2.], 'run2': [3., 4.]})
    save_table_dataset('data.csv', kind, data)
    profiles = []
    for i in range(count):
        profiles.append((f'Profile {i + 1}', pd.DataFrame({
            'sample': ['001', 'run2'], 'condition': ['NA', '=literal'],
            'pair_id': ['007', '007'] if i == 0 else ['', ''],
        })))
    save_uploaded_metadata_profiles(kind, 'input.csv', profiles, active_profile_name=profiles[-1][0])
    return data, profiles


@pytest.mark.parametrize('kind,count', [('protein', 1), ('protein', 2), ('protein', 5), ('phospho', 2)])
def test_all_profiles_round_trip_through_api(kind, count):
    _, profiles = setup_profiles(kind, count)
    if count > 1:
        set_active_profile(list_auto_generated_profiles()[0])
    with TestClient(app) as client:
        download = client.get(f'/api/annotations/current/{kind}/metadata/download?scope=all')
        assert download.status_code == 200, download.text
        assert 'all_profiles.xlsx' in download.headers['content-disposition']
        workbook = load_workbook(io.BytesIO(download.content))
        exported = list(workbook[PROFILE_SHEET].values)
        assert [row[0] for row in exported[1:]] == [name for name, _ in profiles]
        assert len(exported) == count + 1
        assert workbook['Samples']['C3'].data_type == 's'
        assert workbook['Samples']['C3'].value == '=literal'
        workbook.close()
        reset_all_uploaded_metadata()
        uploaded = client.post('/api/annotations/metadata/upload', data={'kind': kind},
            files={'file': ('metadata.xlsx', download.content, 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')})
        assert uploaded.status_code == 200, uploaded.text
    assert uploaded.json()['profileNames'] == [name for name, _ in profiles]
    for name, expected in profiles:
        pd.testing.assert_frame_equal(get_uploaded_metadata(kind, profile_name=name).frame, expected)
    assert len(list_auto_generated_profiles()) == 2**count - count - 1
    assert get_active_profile() == profiles[0][0]


def test_export_prefers_current_annotation_pairing_and_keeps_selected_csv():
    data, profiles = setup_profiles()
    name, metadata = profiles[0]
    metadata = metadata.copy()
    metadata['pair_id'] = ['edited', 'edited']
    save_annotation(kind='protein', source_data=data, log2_data=data, filtered_data=data,
        metadata=metadata, is_log2_transformed=True, metadata_source='manual',
        filter_config=AnnotationFilterConfig(), auto_detected=False, warnings=[],
        metadata_profile_name=name)
    set_active_profile(name)
    records, active = read_metadata_profiles(io.BytesIO(export_metadata_profiles('protein')), 'protein', data)
    assert active == name
    pd.testing.assert_frame_equal(dict(records)[name], metadata)
    with TestClient(app) as client:
        selected = client.get('/api/annotations/current/protein/metadata/download?scope=selected')
        assert selected.status_code == 200
        assert selected.headers['content-type'].startswith('text/csv')
        pd.testing.assert_frame_equal(pd.read_csv(io.BytesIO(selected.content), dtype=str, keep_default_na=False), metadata)


@pytest.mark.parametrize('corruption', ['sample', 'profile', 'active', 'kind'])
def test_invalid_workbook_does_not_change_existing_profiles(corruption):
    _, profiles = setup_profiles()
    workbook = load_workbook(io.BytesIO(export_metadata_profiles('protein')))
    if corruption == 'sample':
        workbook['Samples']['B2'] = 'missing-run'
    elif corruption == 'profile':
        workbook[PROFILE_SHEET]['A3'] = workbook[PROFILE_SHEET]['A2'].value
    elif corruption == 'active':
        workbook[PROFILE_SHEET]['C2'] = 'yes'
    else:
        workbook[PROFILE_SHEET]['B2'] = 'phospho'
    output = io.BytesIO()
    workbook.save(output)
    before_active = get_active_profile()
    with TestClient(app) as client:
        response = client.post('/api/annotations/metadata/upload', data={'kind': 'protein'},
            files={'file': ('metadata.xlsx', output.getvalue(), 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')})
    assert response.status_code == 400
    assert get_active_profile() == before_active
    for name, frame in profiles:
        pd.testing.assert_frame_equal(get_uploaded_metadata('protein', profile_name=name).frame, frame)


def test_ordinary_excel_metadata_upload_remains_supported():
    setup_profiles(count=1)
    source = pd.DataFrame({'sample': ['001', 'run2'], 'condition': ['A', 'B']})
    output = io.BytesIO()
    source.to_excel(output, index=False)
    with TestClient(app) as client:
        response = client.post('/api/annotations/metadata/upload', data={'kind': 'protein'},
            files={'file': ('ordinary.xlsx', output.getvalue(), 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')})
    assert response.status_code == 200, response.text


def test_malformed_workbook_reports_validation_error_and_rewinds():
    stream = io.BytesIO(b"not an Excel workbook")
    with pytest.raises(ValueError, match="valid .xlsx workbook"):
        read_metadata_profiles(stream, "protein", pd.DataFrame())
    assert stream.tell() == 0


def test_profile_workbook_respects_expanded_byte_budget(monkeypatch):
    data, _ = setup_profiles(count=1)
    content = export_metadata_profiles("protein")
    monkeypatch.setenv("COPYLOT_UPLOAD_MAX_FRAME_BYTES", "100")
    with pytest.raises(ValueError, match="above the"):
        read_metadata_profiles(io.BytesIO(content), "protein", data)


def test_profile_workbook_counts_blank_rows_in_cell_budget(monkeypatch):
    data, _ = setup_profiles(count=1)
    workbook = load_workbook(io.BytesIO(export_metadata_profiles("protein")))
    workbook["Samples"]["A100"].number_format = "@"
    output = io.BytesIO()
    workbook.save(output)
    workbook.close()
    monkeypatch.setenv("COPYLOT_UPLOAD_MAX_CELLS", "40")
    with pytest.raises(ValueError, match="cell limit"):
        read_metadata_profiles(io.BytesIO(output.getvalue()), "protein", data)
