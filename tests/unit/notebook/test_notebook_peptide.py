from notebook_replay_support import notebook_code, replay_cells, set_original_path

import io
import json
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

from annotation_support import save_log2_protein_annotation
from app.services.dataset_store import save_peptide_path
from app.services.peptide_metadata_store import save_peptide_metadata
from app.services.notebook_peptide import render_coverage
from app.services.notebook_archive import build_notebook_archive
from app.services.notebook_snapshots import snapshots
from app.services.plot_images import qc_peptide_coverage_plot, _peptide_coverage_summary
from app.services.session_context import get_session_id


@pytest.mark.parametrize('level,counts,total', [('peptide', [1, 2], 2), ('precursor', [2, 3], 3)])
@pytest.mark.parametrize('include_id', [False, True])
def test_coverage_counts_parity_and_independent_notebook(peptide_input, tmp_path, level, counts, total, include_id):
    options = dict(level=level, include_id=include_id, summary=include_id, text=True, text_size=12, random_seed=31, quantity_column='Precursor.Quantity', dpi=80,
                   advanced={'conditionColors': {'A': '#771133'}, 'sampleLabels': {'run1': 'First run'},
                             'conditionLabels': {'B': 'Group B'}}, reference_line=False)
    expected = qc_peptide_coverage_plot(**options)
    actual, token = render_coverage(**options)
    assert actual == expected
    recipe = snapshots.get(get_session_id(), token)
    table = pd.read_parquet(io.BytesIO(recipe.resources['table.parquet']))
    assert table.Value.tolist() == counts
    assert recipe.provenance['totalQuantifiedFeatures'] == total
    assert recipe.provenance['featureColumn'] == ('Stripped.Sequence' if level == 'peptide' else 'Precursor.Id')
    assert not recipe.uses_imputed_input
    original_table, original_total, _ = _peptide_coverage_summary(include_id, 'Precursor.Quantity', level)
    pd.testing.assert_frame_equal(table, original_table)
    assert total == original_total
    _, payload = build_notebook_archive(recipe)
    output = tmp_path / 'export'
    output.mkdir()
    with ZipFile(io.BytesIO(payload)) as archive:
        for name in archive.namelist():
            assert str(peptide_input).encode() not in archive.read(name)
        archive.extractall(output)
        name = next(name for name in archive.namelist() if name.endswith('.ipynb'))
    notebook = json.loads((output / name).read_text())
    cells = notebook_code(notebook)
    cells = set_original_path(cells, output)
    result = replay_cells(cells, output, timeout=60)
    assert result.returncode == 0, result.stderr
    assert (output / 'reproduced.png').read_bytes() == actual
    assert (output / 'reproduced_table.csv').exists()


def test_old_quantity_override_does_not_change_snapshot_input(peptide_input):
    _, first = render_coverage(quantity_column='Precursor.Quantity', dpi=80)
    _, second = render_coverage(quantity_column='Alternative.Quantity', dpi=80)
    first_recipe = snapshots.get(get_session_id(), first)
    second_recipe = snapshots.get(get_session_id(), second)
    assert first_recipe.provenance['totalQuantifiedFeatures'] == 2
    assert second_recipe.provenance['totalQuantifiedFeatures'] == 2
    assert second_recipe.settings['quantity_column'] == 'Precursor.Quantity'
    assert snapshots.get(get_session_id(), first).resources == first_recipe.resources


@pytest.fixture
def rt_input(tmp_path):
    frame = pd.DataFrame({'Predicted.RT': ['10', '12', '13', 'bad', '16', '20'],
        'RT': [11., 13., 15., 14., np.nan, 18.],
        'File.Name': [r'C:\notebook_private_rt_origin\run.raw'] * 6})
    path = tmp_path / 'rt.csv'
    frame.to_csv(path, index=False)
    save_peptide_path(str(path))
    return frame


@pytest.mark.parametrize('method', ['Scatter Plot', 'Density Plot', 'Hexbin Plot'])
@pytest.mark.parametrize('styled', [False, True])
def test_rt_exact_image_and_isolated_replay(rt_input, tmp_path, method, styled):
    from app.services.notebook_peptide import render_rt
    from app.services.peptide_tools import peptide_rt_plot
    options = dict(method=method, bins=20, dpi=80, add_line=styled, header=not styled,
                   sampling_seed=21, advanced={'xTickLabelSize': 8, 'yTickLabelSize': 9} if styled else None)
    expected = peptide_rt_plot(**options)
    image, token = render_rt(**options)
    assert image == expected
    recipe = snapshots.get(get_session_id(), token)
    assert not recipe.uses_imputed_input
    captured = pd.read_parquet(io.BytesIO(recipe.resources['input.parquet']))
    pd.testing.assert_frame_equal(captured, rt_input[['Predicted.RT', 'RT']])
    _, blob = build_notebook_archive(recipe)
    with ZipFile(io.BytesIO(blob)) as archive:
        assert all(b'notebook_private_rt_origin' not in archive.read(name) for name in archive.namelist())
        archive.extractall(tmp_path)
        name = next(name for name in archive.namelist() if name.endswith('.ipynb'))
    notebook = json.loads((tmp_path / name).read_text())
    cells = notebook_code(notebook)
    cells = set_original_path(cells, tmp_path)
    replay = replay_cells(cells, tmp_path, timeout=60)
    assert replay.returncode == 0, replay.stderr
    assert (tmp_path / 'reproduced.png').read_bytes() == image


def test_rt_large_scatter_repeats_sampling(rt_input, monkeypatch):
    from app.services.notebook_peptide import render_rt
    from app.services import peptide_tools
    rng = np.random.default_rng(82)
    data = pd.DataFrame(rng.normal(30, 10, (100001, 2)), columns=['Predicted.RT', 'RT'])
    monkeypatch.setattr(peptide_tools, 'get_peptide_frame', lambda **kwargs: data.copy())
    options = dict(method='Scatter Plot', dpi=72, sampling_seed=17, add_line=True)
    image, token = render_rt(**options)
    assert image == peptide_tools.peptide_rt_plot(**options)
    assert snapshots.get(get_session_id(), token).settings['sampling_seed'] == 17


def test_rt_missing_columns_and_empty_pairs_keep_original_errors(rt_input, monkeypatch):
    from app.services.notebook_peptide import render_rt
    from app.services import peptide_tools
    for frame, message in [(pd.DataFrame({'RT': [1]}), 'Required peptide column'),
                           (pd.DataFrame({'RT': [1], 'Predicted.RT': ['bad']}), 'No numeric')]:
        monkeypatch.setattr(peptide_tools, 'get_peptide_frame', lambda **kwargs: frame.copy())
        with pytest.raises(ValueError, match=message):
            render_rt()


def test_rt_route_snapshot_is_immutable_and_owned(rt_input, monkeypatch):
    from fastapi import FastAPI
    from fastapi.testclient import TestClient
    from app.api.routes_plots import router
    from app.services.session_context import set_session_id, reset_session_id
    monkeypatch.setenv('COPYLOT_SESSION_ISOLATION', '1')
    app = FastAPI()
    app.include_router(router)
    with TestClient(app) as client:
        response = client.get('/api/plots/peptide/rt.png', params={'bins': 20, 'dpi': 80})
        assert response.status_code == 200
        url = '/api/plots/notebooks/' + response.headers['X-CoPYlot-Notebook'] + '.zip'
        first = client.get(url)
        assert first.status_code == 200
        assert client.get(url).content == first.content
        other = set_session_id('other-rt-owner')
        try:
            assert client.get(url).status_code == 404
        finally:
            reset_session_id(other)


@pytest.fixture
def modification_input(tmp_path):
    from app.services.dataset_store import reset_all_datasets
    from app.services.annotation_store import reset_all_annotations
    reset_all_datasets()
    reset_all_annotations()
    names = [r'C:\modification_private_origin\ctrl_01.raw', r'C:\modification_private_origin\treated_01.raw']
    sequences = ['AC(UniMod:4)M(UniMod:35)', 'AC(UniMod:4)M(UniMod:35)',
        '(UniMod:1)PEPTIDE', 'M(UniMod:35)M(UniMod:35)', 'UNMODIFIED', 'M(UniMod:350)']
    frame = pd.DataFrame({'File.Name': [names[0]] * 6 + [names[1]] * 6,
        'Modified.Sequence': sequences * 2,
        'Precursor.Quantity': [1., 3., 0., np.inf, 7., 10., 0., 0., 4., 2., 9., 10.]})
    path = tmp_path / 'modifications.csv'
    frame.to_csv(path, index=False)
    save_peptide_path(str(path))
    metadata = pd.DataFrame({'sample': names, 'condition': ['A', 'B']})
    save_peptide_metadata('metadata.csv', metadata)
    yield metadata
    reset_all_datasets()
    reset_all_annotations()


@pytest.mark.parametrize('protein_metadata', [False, True])
@pytest.mark.parametrize('include_id', [False, True])
def test_modification_counts_and_independent_plot_replay(modification_input, tmp_path, protein_metadata, include_id):
    from app.services.notebook_peptide import render_modification
    from app.services.peptide_tools import peptide_modification_plot
    if protein_metadata:
        from app.services.dataset_store import save_table_dataset
        data = pd.DataFrame({'ProteinNames': ['P1'], **{name: [10.] for name in modification_input['sample']}})
        save_table_dataset('protein.csv', 'protein', data)
        save_log2_protein_annotation(data, modification_input)
    options = dict(include_id=include_id, dpi=80, header=not include_id, legend=not include_id,
                   advanced={'xTickLabelSize': 9, 'yTickLabelSize': 8})
    image, token = render_modification(**options)
    assert image == peptide_modification_plot(**options)
    recipe = snapshots.get(get_session_id(), token)
    assert recipe.provenance['metadataSource'] == ('protein annotation' if protein_metadata else 'peptide annotation')
    counts = pd.read_parquet(io.BytesIO(recipe.resources['counts.parquet']))
    assert counts['Carbamidomethylation'].tolist() == [1, 0]
    assert counts['Oxidation'].tolist() == [1, 1]
    assert counts['Acetylation'].tolist() == [0, 1]
    assert counts['UniqueTrackedModifiedPeptides'].tolist() == [1, 2]
    _, blob = build_notebook_archive(recipe)
    with ZipFile(io.BytesIO(blob)) as archive:
        assert all(b'modification_private_origin' not in archive.read(name) for name in archive.namelist())
        archive.extractall(tmp_path)
        name = next(name for name in archive.namelist() if name.endswith('.ipynb'))
    notebook = json.loads((tmp_path / name).read_text())
    cells = notebook_code(notebook)
    cells = set_original_path(cells, tmp_path)
    replay = replay_cells(cells, tmp_path, timeout=60)
    assert replay.returncode == 0, replay.stderr
    assert (tmp_path / 'reproduced.png').read_bytes() == image
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / 'reproduced_modification_counts.csv', index_col=0), counts)


def test_modification_no_detected_classes_remains_exportable(modification_input, monkeypatch):
    from app.services import peptide_tools
    from app.services.notebook_peptide import render_modification
    frame = peptide_tools.get_peptide_frame()
    frame['Precursor.Quantity'] = 0.
    monkeypatch.setattr(peptide_tools, 'get_peptide_frame', lambda **kwargs: frame.copy())
    image, token = render_modification(dpi=80)
    assert image == peptide_tools.peptide_modification_plot(dpi=80)
    recipe = snapshots.get(get_session_id(), token)
    assert pd.read_parquet(io.BytesIO(recipe.resources['counts.parquet'])).to_numpy().sum() == 0


def test_modification_sample_mismatch_keeps_original_error(modification_input):
    from app.services import peptide_tools
    from app.services.notebook_peptide import render_modification
    save_peptide_metadata('unmatched.csv', pd.DataFrame({'sample': ['other_99'], 'condition': ['C']}))
    # Both uploaded runs end in 01; fallback must reject that ambiguous match.
    for renderer in (peptide_tools.peptide_modification_plot, render_modification):
        with pytest.raises(ValueError, match='ambiguous'):
            renderer(dpi=80)


@pytest.fixture
def cleavage_input(modification_input, monkeypatch):
    from app.services import peptide_tools
    sequences = ['AKPAK', 'AKPAK', 'AKRAK', 'PEPTIDE', 'RPAK', 'AK']
    frame = pd.DataFrame({'File.Name': [modification_input['sample'][0]] * 6 + [modification_input['sample'][1]] * 6,
        'Stripped.Sequence': sequences * 2,
        'Precursor.Quantity': [1., 2., 3., 4., 0., np.inf, 5., 0., 0., 2., 3., 0.]})
    monkeypatch.setattr(peptide_tools, 'get_peptide_frame', lambda **kwargs: frame.copy())
    return frame


@pytest.mark.parametrize('protease', ['trypsin', 'trypsin_lysc'])
@pytest.mark.parametrize('styled', [False, True])
def test_cleavage_exact_counts_and_independent_replay(cleavage_input, tmp_path, protease, styled):
    from app.services.notebook_peptide import render_missed_cleavage
    from app.services.peptide_tools import peptide_missed_cleavage_plot
    options = dict(protease=protease, dpi=80, include_id=styled, text=not styled, header=not styled,
        text_size=11, advanced={'xTickLabelSize': 8, 'yTickLabelSize': 9} if styled else None)
    image, token = render_missed_cleavage(**options)
    assert image == peptide_missed_cleavage_plot(**options)
    recipe = snapshots.get(get_session_id(), token)
    summary = pd.read_parquet(io.BytesIO(recipe.resources['summary.parquet']))
    assert summary.Total.eq(3).all()
    np.testing.assert_allclose(summary.groupby('Sample').Percentage.sum(), [100., 100.])
    a = summary[summary.Sample.str.startswith('A')].set_index('Count').Occurrences.to_dict()
    b = summary[summary.Sample.str.startswith('B')].set_index('Count').Occurrences.to_dict()
    assert a == ({0: 2, 2: 1} if protease == 'trypsin' else {0: 1, 1: 1, 2: 1})
    assert b == ({0: 3} if protease == 'trypsin' else {0: 2, 1: 1})
    _, blob = build_notebook_archive(recipe)
    with ZipFile(io.BytesIO(blob)) as archive:
        assert all(b'modification_private_origin' not in archive.read(name) for name in archive.namelist())
        archive.extractall(tmp_path)
        name = next(name for name in archive.namelist() if name.endswith('.ipynb'))
    notebook = json.loads((tmp_path / name).read_text())
    cells = notebook_code(notebook)
    cells = set_original_path(cells, tmp_path)
    replay = replay_cells(cells, tmp_path, timeout=60)
    assert replay.returncode == 0, replay.stderr
    assert (tmp_path / 'reproduced.png').read_bytes() == image
    pd.testing.assert_frame_equal(pd.read_csv(tmp_path / 'reproduced_missed_cleavages.csv'), summary)


def test_cleavage_empty_and_invalid_protease_preserve_errors(cleavage_input):
    from app.services.notebook_peptide import render_missed_cleavage
    from app.services.peptide_tools import peptide_missed_cleavage_plot
    for renderer in (peptide_missed_cleavage_plot, render_missed_cleavage):
        with pytest.raises(ValueError, match='Protease must'):
            renderer(protease='unknown')
    cleavage_input['Precursor.Quantity'] = 0.
    for renderer in (peptide_missed_cleavage_plot, render_missed_cleavage):
        with pytest.raises(ValueError, match='No peptide intensities'):
            renderer()
