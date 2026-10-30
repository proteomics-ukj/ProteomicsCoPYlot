"""Preparation reuse must not change results, ownership or invalidation."""
import os
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from app.services import peptide_preparation_cache as cache
from app.services import peptide_tools as peptides
from app.services.dataset_store import save_peptide_path
from app.services.peptide_metadata_store import save_peptide_metadata
from app.services.runtime_cache import invalidate_runtime_cache


@pytest.fixture
def source(tmp_path):
    path = tmp_path / 'peptides.csv'
    pd.DataFrame({'File.Name': ['run1', 'run2'], 'Modified.Sequence': ['AC(UniMod:4)', 'PEPTIDE'],
                  'Stripped.Sequence': ['ACK', 'PEPTIDE'], 'Precursor.Quantity': [10., 20.],
                  'ProteinNames': ['P1', 'P2'], 'Predicted.RT': [1., 2.], 'RT': [2., 3.]}).to_csv(path, index=False)
    save_peptide_path(str(path))
    save_peptide_metadata('metadata.csv', pd.DataFrame({'sample': ['run1', 'run2'], 'condition': ['A', 'B']}))
    return path


def test_mutation_revision_file_and_session_isolation(source, monkeypatch):
    factory = Mock(side_effect=lambda: {'frame': pd.DataFrame({'x': [1]})})
    first = cache.cached_peptide_preparation('test', (), factory)
    first['frame'].iloc[0, 0] = 999
    assert cache.cached_peptide_preparation('test', (), factory)['frame'].iloc[0, 0] == 1
    assert factory.call_count == 1
    invalidate_runtime_cache('changed-profile-or-filter')
    cache.cached_peptide_preparation('test', (), factory)
    assert factory.call_count == 2
    stamp = source.stat()
    os.utime(source, ns=(stamp.st_atime_ns, stamp.st_mtime_ns + 1_000_000))
    cache.cached_peptide_preparation('test', (), factory)
    assert factory.call_count == 3
    monkeypatch.setattr(cache, 'get_session_id', lambda: 'another-session')
    cache.cached_peptide_preparation('test', (), factory)
    assert factory.call_count == 4


def test_cache_budget_evicts_and_oversized_results_still_work(source, monkeypatch):
    monkeypatch.setattr(cache, '_CACHE', cache.OrderedDict())
    monkeypatch.setattr(cache, '_MAX_BYTES', 180)
    one = Mock(return_value=b'x' * 100)
    two = Mock(return_value=b'y' * 100)
    cache.cached_peptide_preparation('one', (), one)
    cache.cached_peptide_preparation('two', (), two)
    cache.cached_peptide_preparation('one', (), one)
    assert one.call_count == 2
    big = Mock(return_value=b'z' * 1000)
    assert cache.cached_peptide_preparation('big', (), big) == b'z' * 1000
    cache.cached_peptide_preparation('big', (), big)
    assert big.call_count == 2
    assert sum(entry[1] for entry in cache._CACHE.values()) <= 180


@pytest.mark.parametrize('kind', ['modification', 'missed_cleavage'])
def test_restyle_reuses_preparation_and_notebook_inputs(source, monkeypatch, kind):
    from app.services import notebook_peptide as notebooks
    from app.services.notebook_snapshots import snapshots
    from app.services.session_context import get_session_id
    module = getattr(notebooks, f'_{kind}_module')()
    prepare = Mock(wraps=module.prepare_input)
    monkeypatch.setattr(module, 'prepare_input', prepare)
    summary_name = 'peptide_modification_counts' if kind == 'modification' else 'peptide_missed_cleavage_summary'
    summary = Mock(wraps=getattr(module, summary_name))
    monkeypatch.setattr(module, summary_name, summary)
    render = getattr(notebooks, f'render_{kind}')
    first, first_token = render(dpi=72)
    restyled, second_token = render(dpi=72, header=False)
    repeated, third_token = render(dpi=72)
    assert first == repeated and first != restyled
    assert (prepare.call_count, summary.call_count) == (1, 1)
    assert all([first_token, second_token, third_token])
    recipe = snapshots.get(get_session_id(), third_token)
    assert recipe.resources['original.png'] == first
    assert 'input.parquet' in recipe.resources
    save_peptide_metadata('changed.csv', pd.DataFrame({'sample': ['run1', 'run2'], 'condition': ['C', 'D']}))
    render(dpi=72)
    assert prepare.call_count == summary.call_count == 2


def test_index_matches_original_group_order_and_detection(monkeypatch):
    frame = pd.DataFrame({'ProteinNames': [' P1 ;P2', 'P2', 'P1', None, 'P1', 'P2', 'P2'],
                          'Stripped.Sequence': ['ABC', 'DEF', 'ABC', 'GHI', 'ZERO', None, 'INF'],
                          'Precursor.Quantity': [1, 2, 3, 4, 0, 5, np.inf]})
    monkeypatch.setattr(peptides, 'get_peptide_frame', lambda: frame)
    index = peptides._protein_peptide_index()
    for protein in ['P1', 'P2', '', 'missing']:
        assert index.get(protein, []) == peptides._protein_peptides(frame, protein)


def test_column_projection_does_not_mutate_shared_input(source):
    selected = peptides.get_peptide_frame(columns=['RT'])
    selected.iloc[0, 0] = 100
    assert list(selected.columns) == ['RT']
    assert peptides.get_peptide_frame()['RT'].tolist() == [2., 3.]


def test_external_rewrite_and_output_formats_do_not_reuse_old_results(source):
    from app.services.notebook_peptide import _cached_image
    from app.services.runtime_cache import runtime_cache_version
    assert peptides.get_peptide_frame()['RT'].iloc[0] == 2.
    frame = pd.read_csv(source)
    frame['RT'] = [9., 10.]
    frame.to_csv(source, index=False)
    assert peptides.get_peptide_frame()['RT'].tolist() == [9., 10.]
    revision = runtime_cache_version()
    assert _cached_image('format-test', {}, 'png', revision, lambda: b'png') == b'png'
    assert _cached_image('format-test', {}, 'svg', revision, lambda: b'svg') == b'svg'
    assert _cached_image('format-test', {}, 'png', revision + 1, lambda: b'changed') == b'changed'
