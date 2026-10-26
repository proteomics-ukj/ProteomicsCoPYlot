import io
import json
import gzip
from zipfile import ZipFile

import pandas as pd

from app.api.routes_annotations import generate_annotation, normalize_annotation
from app.schemas.annotation import AnnotationGenerateRequest, AnnotationNormalizeRequest
from app.services.annotation_store import get_annotation
from app.services.dataset_store import save_table_dataset, save_peptide_path


def test_protein_phospho_and_normalization_keep_independent_filter_settings():
    frame = pd.DataFrame({"ProteinNames": ["ALBU_HUMAN", "AMY1A_HUMAN", "SAMPLE", "ALBU_HUMAN;SAMPLE"],
                          "A": [12., 15., 18., 22.], "B": [13., 14., 19., 23.]})
    conditions = [{"name": "one", "columns": ["A", "B"]}]
    for kind in ['protein', 'phospho']:
        save_table_dataset('input.csv', kind, frame)
        response = generate_annotation(AnnotationGenerateRequest(kind=kind, isLog2Transformed=True,
            conditions=conditions, filter={"contaminants": {"enabled": kind == 'protein', "categories": [4]}}))
        assert response.filteredRows == (3 if kind == 'protein' else 4)
    assert get_annotation('protein').filtered_data.ProteinNames.tolist() == ['AMY1A_HUMAN', 'SAMPLE', 'ALBU_HUMAN;SAMPLE']
    response = normalize_annotation(AnnotationNormalizeRequest(kind='protein', normalization='median'))
    assert response.filteredRows == 3
    assert response.filter.contaminants.categories == [4]
    assert any('Metadata SHA-256' in w for w in response.warnings)
    assert not get_annotation('phospho').filter_config.contaminants.enabled


def test_peptide_export_replays_filter_from_original_upload(tmp_path):
    from app.services.notebook_peptide import render_rt
    from app.services.notebook_snapshots import snapshots
    from app.services.session_context import get_session_id
    from app.services.notebook_archive import build_notebook_archive
    protein = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN', 'SAMPLE'], 'A': [10., 20.]})
    save_table_dataset('protein.csv', 'protein', protein)
    generate_annotation(AnnotationGenerateRequest(kind='protein', isLog2Transformed=True,
        conditions=[{'name': 'one', 'columns': ['A']}], filter={'contaminants': {'enabled': True}}))
    peptides = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN', 'SAMPLE', 'ALBU_HUMAN;SAMPLE', 'SAMPLE2'],
        'Predicted.RT': [10., 20., 30., 40.], 'RT': [11., 21., 31., 41.], 'File.Name': ['run'] * 4})
    path = tmp_path / 'peptides.csv'
    peptides.to_csv(path, index=False)
    save_peptide_path(str(path))
    _, token = render_rt(method='Scatter Plot', dpi=80)
    recipe = snapshots.get(get_session_id(), token)
    config = json.loads(recipe.resources['processing.json'])
    assert len(config['contaminant_reference']['records']) == 125
    assert config['contaminants']['enabled']
    namespace = {}
    exec(recipe.source.code, namespace)
    replayed = namespace['filter_peptide_contaminants'](peptides, config['contaminants'], config['contaminant_reference'])
    assert replayed.ProteinNames.tolist() == ['SAMPLE', 'ALBU_HUMAN;SAMPLE', 'SAMPLE2']
    assert 'filter_peptide_contaminants' in recipe.preparation
    # The downloadable archive includes the exact lookup used for this render.
    _, payload = build_notebook_archive(recipe)
    with ZipFile(io.BytesIO(payload)) as archive:
        assert any(name.endswith('processing.json') for name in archive.namelist())


def test_contaminant_settings_and_rows_survive_session_restore(tmp_path):
    from app.services import restore_state
    from app.services.annotation_store import reset_all_annotations
    from app.services.dataset_store import reset_all_datasets
    frame = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN', 'AMY1A_HUMAN', 'SAMPLE'], 'A': [12., 14., 18.]})
    path = tmp_path / 'protein.csv'
    frame.to_csv(path, index=False)
    save_table_dataset(str(path), 'protein', frame, source_path=str(path))
    generate_annotation(AnnotationGenerateRequest(kind='protein', isLog2Transformed=True,
        conditions=[{'name': 'one', 'columns': ['A']}],
        filter={'contaminants': {'enabled': True, 'categories': [4], 'matchAccessions': False}}))
    before = get_annotation('protein').filtered_data.copy()
    payload = json.loads(gzip.decompress(restore_state.export_state_blob()))
    reset_all_annotations()
    reset_all_datasets()
    restore_state.apply_restore_state(payload, {})
    restored = get_annotation('protein')
    assert restored.filter_config.contaminants.categories == [4]
    assert not restored.filter_config.contaminants.matchAccessions
    pd.testing.assert_frame_equal(restored.filtered_data, before)
