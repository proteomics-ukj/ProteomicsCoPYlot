import pandas as pd
import pytest
from pydantic import ValidationError

from app.schemas.annotation import AnnotationFilterConfig, ContaminantFilterConfig
from app.services.contaminant_filter import contaminant_mask, contaminant_reference
from app.services.annotation_processor import filter_then_normalize


def test_exact_groups_categories_and_missing_identifiers():
    frame = pd.DataFrame({"ProteinNames": [
        "ALBU_HUMAN", "Albumin", "ALBU_HUMAN;ADH1_YEAST", "ALBU_HUMAN;MY_SAMPLE",
        "ALB", "prefix_ALBU_HUMAN", None, "ALBU_HUMAN;", "AMY1A_HUMAN", "ENDOP_GLUC",
        "sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens", "P02768-2", "000000",
    ]})
    reference = contaminant_reference()
    assert not contaminant_mask(frame).any()
    assert contaminant_mask(frame, {"enabled": True}, reference).tolist() == [
        True, True, True, False, False, False, False, False, True, True, True, False, False,
    ]
    selected = contaminant_mask(frame, {"enabled": True, "categories": [4]}, reference)
    assert selected.iloc[:4].tolist() == [True, True, False, False]
    union = contaminant_mask(frame, {"enabled": True, "categories": [3, 4]}, reference)
    assert union.iloc[2]


def test_accession_fallback_is_exact_optional_and_group_complete():
    frame = pd.DataFrame({"ProteinNames": [None, "display name", "ALBU_HUMAN;MY_SAMPLE", "ALBU_HUMAN", "ALBU_HUMAN"],
                          "ProteinAccessions": ["P02768", "P02768", "P02768", "P02768;UNKNOWN", "P02768-2"]})
    assert contaminant_mask(frame, {"enabled": True}, contaminant_reference()).tolist() == [True, True, False, False, False]
    assert contaminant_mask(frame, {"enabled": True, "matchAccessions": False}, contaminant_reference()).tolist() == [False, False, False, True, True]
    assert contaminant_mask(pd.DataFrame({"Protein_group": ["P02768", "P02768;unknown"]}),
                            {"enabled": True}, contaminant_reference()).tolist() == [True, False]
    with pytest.raises(ValueError, match="Gene names are not used"):
        contaminant_mask(pd.DataFrame({"GeneNames": ["ALB"]}), {"enabled": True}, contaminant_reference())
    commercial = pd.DataFrame({'ProteinNames': ['ENDOP_GLUC', 'SAMPLE', 'ENDOP_GLUC;SAMPLE'],
                               'ProteinAccessions': ['000000', '000000', '000000']})
    assert contaminant_mask(commercial, {'enabled': True}, contaminant_reference()).tolist() == [True, False, False]
    commercial['ProteinAccessions'] = 'CON__000000'
    assert contaminant_mask(commercial, {'enabled': True}, contaminant_reference()).tolist() == [True, False, False]


def test_search_engine_prefixes_require_database_matches_and_keep_decoys():
    from app.services.dataset_reader import rename_cols
    from app.services.notebook_preprocessing import preprocessing_source
    labels = ['CON__P02768', 'contam_sp|P02768|ALBU_HUMAN', 'CON__UNKNOWN',
              'CON__P02768;SAMPLE', 'REV__P02768', 'rev_sp|P02768|ALBU_HUMAN',
              'CON__P02768;contam_sp|P02768|ALBU_HUMAN', 'SAMPLE', 'contam_REV__P02768']
    frame = rename_cols(pd.DataFrame({'Protein IDs': labels, 'Majority protein IDs': labels,
                                     'Potential contaminant': ['+'] * len(labels)}))
    reference = contaminant_reference()
    config = {'enabled': True, 'categories': [4]}
    expected = [True, True, False, False, False, False, True, False, False]
    assert contaminant_mask(frame, config, reference).tolist() == expected
    assert not contaminant_mask(frame, {'enabled': True, 'categories': [2]}, reference).any()
    assert contaminant_mask(frame, {**config, 'matchAccessions': False}, reference).tolist() == [
        False, True, False, False, False, False, False, False, False]
    namespace = {}
    exec(preprocessing_source().code, namespace)
    assert namespace['contaminant_mask'](frame, config, reference).tolist() == expected


def test_multi_category_record_matches_any_selected_category():
    from app.services.notebook_preprocessing import preprocessing_source
    namespace = {}
    exec(preprocessing_source().code, namespace)
    reference = contaminant_reference()
    record = next(r for r in reference['records'] if r['entry_name'] == 'ALBU_HUMAN')
    record['gpm_category_id'] = '1;4'
    frame = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN', 'ALBU_HUMAN;ADH1_YEAST']})
    for categories in ([1], [4], [1, 4]):
        assert contaminant_mask(frame, {'enabled': True, 'categories': categories}, reference).tolist() == [True, False]
        assert namespace['contaminant_mask'](frame, {'enabled': True, 'categories': categories}, reference).tolist() == [True, False]
    assert not contaminant_mask(frame, {'enabled': True, 'categories': [2]}, reference).any()


def test_old_sessions_default_off_and_empty_selection_rejected():
    assert not AnnotationFilterConfig().contaminants.enabled
    with pytest.raises(ValidationError):
        ContaminantFilterConfig(enabled=True, categories=[])
    with pytest.raises(ValidationError):
        ContaminantFilterConfig(enabled=True, categories=[6])


def test_original_cambridge_headers_match_without_accessions():
    from app.services.reference_db import resolve_reference_db
    headers = [line for line in resolve_reference_db('crap_ccp.fasta', 'Cambridge FASTA').read_text().splitlines()
               if line.startswith('>')]
    assert len(headers) == 125
    assert contaminant_mask(pd.DataFrame({'ProteinNames': headers}),
        {'enabled': True, 'matchAccessions': False}, contaminant_reference()).all()


def test_contaminant_only_filter_is_described_in_reports():
    from types import SimpleNamespace
    from app.services.report_function import _filter_description, _filter_not_applied_description
    annotation = SimpleNamespace(filter_config=AnnotationFilterConfig(contaminants={"enabled": True}),
        source_data=pd.DataFrame({'A': [1, 2]}), filtered_data=pd.DataFrame({'A': [2]}))
    assert 'Cambridge cRAP contaminant filtering' in _filter_description(annotation)
    assert 'every group member must match' in _filter_description(annotation)
    assert '10.5281/zenodo.15115102' in _filter_description(annotation)
    assert 'Cambridge Centre for Proteomics / camprotR' in _filter_description(annotation)
    assert 'Not applied to this plot' in _filter_not_applied_description(annotation)


def test_all_rows_removed_can_be_filtered_but_not_normalized():
    frame = pd.DataFrame({'ProteinNames': ['ALBU_HUMAN'], 'A': [12.]})
    metadata = pd.DataFrame({'sample': ['A'], 'condition': ['one']})
    args = dict(min_present=0, filter_mode='overall', contaminants={'enabled': True})
    assert filter_then_normalize(frame, metadata, True, 'none', **args)[1].empty
    with pytest.raises(ValueError, match='no rows remain'):
        filter_then_normalize(frame, metadata, True, 'median', **args)


@pytest.mark.parametrize("fit", ["filtered_rows", "complete_features"])
def test_order_normalization_and_portable_replay(fit):
    from app.services.notebook_preprocessing import preprocessing_source
    raw = pd.DataFrame({"ProteinNames": ["ALBU_HUMAN", "AMY1A_HUMAN", "SAMPLE", "SAMPLE2"],
                        "A": [100., None, 10., 20.], "B": [200., 12., 11., 22.]})
    metadata = pd.DataFrame({"sample": ["A", "B"], "condition": ["one", "one"]})
    filtering = {"minPresent": 2, "mode": "overall", "method": "at_least_n", "minPercentage": 70,
                 "normalizationFit": fit, "contaminants": {"enabled": True}}
    reference = contaminant_reference()
    all_rows, filtered = filter_then_normalize(raw, metadata, True, "median", min_present=2,
        filter_mode="overall", normalization_fit=fit, contaminants={"enabled": True}, contaminant_snapshot=reference)
    assert len(all_rows) == 4
    assert filtered.ProteinNames.tolist() == ["SAMPLE", "SAMPLE2"]
    assert "removed 1 rows after" in filtered.attrs['contaminant_filter']
    assert filtered.A.median() == pytest.approx(filtered.B.median())
    namespace = {}
    exec(preprocessing_source().code, namespace)
    reproduced = namespace['replay_preprocessing'](raw, metadata, {
        "source": "filtered", "kind": "protein", "is_log2": True, "normalization": "median",
        "filter": filtering, "contaminant_reference": reference, "linear_output": False})
    pd.testing.assert_frame_equal(reproduced, filtered)


def test_peptides_inherit_protein_only_and_preserve_original(monkeypatch, tmp_path):
    from types import SimpleNamespace
    from app.services import peptide_tools
    raw = pd.DataFrame({"ProteinNames": ["ALBU_HUMAN", "ALBU_HUMAN;SAMPLE", "SAMPLE"], "Quantity": [1., 2., 3.]})
    path = tmp_path / 'peptides.tsv'
    raw.to_csv(path, sep='\t', index=False)
    monkeypatch.setattr(peptide_tools, 'get_current_dataset', lambda kind: SimpleNamespace(path=path))
    calls = []
    def annotation(kind):
        calls.append(kind)
        return SimpleNamespace(filter_config=AnnotationFilterConfig(contaminants={"enabled": True}))
    monkeypatch.setattr(peptide_tools, 'get_annotation', annotation)
    assert peptide_tools.get_peptide_frame().ProteinNames.tolist() == ["ALBU_HUMAN;SAMPLE", "SAMPLE"]
    assert calls == ['protein']
    assert len(peptide_tools.get_peptide_frame(apply_contaminants=False)) == 3
