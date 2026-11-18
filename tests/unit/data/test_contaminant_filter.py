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
    assert contaminant_mask(commercial, {'enabled': True, 'categories': [1, 2, 3, 4, 5]}, contaminant_reference()).tolist() == [True, False, False]


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


@pytest.mark.parametrize("config, removed", [
    ({"enabled": False}, []),
    ({"enabled": True, "categories": [4]}, ["CON__P02768", "P02768"]),
    ({"enabled": True, "categories": [2]}, []),
    ({"enabled": True, "categories": [4], "matchAccessions": False}, []),
    ({"enabled": True, "categories": [6]}, ["CON__P02768", "CON__UNKNOWN", "CON__P02768;SAMPLE"]),
    ({"enabled": True, "categories": [4, 6]}, ["CON__P02768", "CON__UNKNOWN", "P02768", "CON__P02768;SAMPLE"]),
])
def test_maxquant_upload_flags_and_database_filter_have_distinct_effects(config, removed):
    import io
    from app.services.dataset_reader import read_dataframe

    names = ["CON__P02768", "CON__UNKNOWN", "P02768", "CON__P02768;SAMPLE", "SITE_ONLY", "REV__P02768"]
    raw = pd.DataFrame({"Protein IDs": names, "Majority protein IDs": names,
                        "Potential contaminant": ["+", "+", "", "+", "", ""],
                        "Reverse": ["", "", "", "", "", "+"],
                        "Only identified by site": ["", "", "", "", "+", ""],
                        "LFQ intensity a": [100.] * 6, "LFQ intensity b": [200.] * 6})
    frame = read_dataframe("proteinGroups.txt", io.BytesIO(raw.to_csv(sep="\t", index=False).encode()))
    # Decoys and groups identified only by a modification site are removed on loading.
    loaded = names[:4]
    assert frame.ProteinNames.tolist() == loaded
    assert frame["Potential contaminant"].fillna("").tolist() == ["+", "+", "", "+"]
    metadata = pd.DataFrame({"sample": ["LFQ intensity a", "LFQ intensity b"], "condition": ["one", "one"]})
    original, filtered = filter_then_normalize(frame, metadata, False, "none", min_present=0,
                                               filter_mode="overall", contaminants=config)
    assert original.ProteinNames.tolist() == loaded
    assert filtered.ProteinNames.tolist() == [name for name in loaded if name not in removed]


def test_old_sessions_default_off_and_empty_selection_rejected():
    assert not AnnotationFilterConfig().contaminants.enabled
    with pytest.raises(ValidationError):
        ContaminantFilterConfig(enabled=True, categories=[])
    with pytest.raises(ValidationError):
        ContaminantFilterConfig(enabled=True, categories=[7])


@pytest.mark.parametrize("categories", [[6], [4, 6]])
@pytest.mark.parametrize("accessions", [False, True])
def test_upstream_category_is_independent_and_replays_portably(categories, accessions):
    from app.services.notebook_preprocessing import preprocessing_source
    from proteomicscopylot.notebook_functions import filter_contaminants
    frame = pd.DataFrame({
        "ProteinNames": ["UNKNOWN", "CON__UNKNOWN", "contam_OTHER", "CON__X;SAMPLE", "CON__X;CON__Y",
                         "REV__X", "SITE_ONLY", "ALBU_HUMAN", "CON__X;", None, "CON__X;SAMPLE"],
        "Potential contaminant": [" + ", "", None, "", "", "", "", "", "", None, "+"],
        "Only identified by site": [""] * 6 + ["+"] + [""] * 4,
    }, index=range(10, 21))
    config = {"enabled": True, "categories": categories, "matchAccessions": accessions}
    expected = [True, True, True, False, True, False, False, 4 in categories, False, False, True]
    reference = contaminant_reference()
    assert contaminant_mask(frame, config, reference).tolist() == expected
    assert not contaminant_mask(frame, {**config, "enabled": False}, reference).any()
    namespace = {}
    exec(preprocessing_source().code, namespace)
    assert namespace["contaminant_mask"](frame, config, reference).tolist() == expected
    pd.testing.assert_frame_equal(filter_contaminants(frame, enabled=True, categories=categories,
                                                     match_accessions=accessions), frame.loc[[not v for v in expected]])


def test_upstream_only_needs_no_database_or_identifier_column():
    frame = pd.DataFrame({"Potential contaminant": ["+", "", None], "A": [1., 2., 3.]})
    assert contaminant_mask(frame, {"enabled": True, "categories": [6]}).tolist() == [True, False, False]


@pytest.mark.parametrize("name_column,accession_column", [("Protein", "Protein ID"), ("PG.ProteinNames", "PG.ProteinAccessions"),
                                                        ("Protein.Names", "Protein.Group")])
def test_pipeline_contaminant_prefixes_survive_identifier_mapping(name_column, accession_column):
    from app.services.dataset_reader import rename_cols
    frame = rename_cols(pd.DataFrame({
        name_column: ["contam_sp|X1|CONTAM_PROTEIN", "CON__X2", "contam_X3;SAMPLE", "CON__X4"],
        accession_column: ["X1", "X2", "X3;SAMPLE", "OTHER"],
    }))
    assert contaminant_mask(frame, {"enabled": True, "categories": [6], "matchAccessions": False}).tolist() == [True, True, False, False]


def test_crap_fasta_prefix_matches_the_database_and_marks_upstream_contaminants():
    """DIA-NN searched with a cRAP FASTA writes cRAP-TRYP_PIG / cRAP-P00761 (PXD081067)."""
    from app.services.dataset_reader import rename_cols
    from app.services.notebook_preprocessing import preprocessing_source
    from proteomicscopylot.notebook_functions import filter_contaminants
    frame = rename_cols(pd.DataFrame({
        "Protein.Names": ["cRAP-ALBU_BOVIN", "K2C1_HUMAN;cRAP-K2C1_HUMAN", "cRAP-UNKNOWN",
                          "cRAP-ALBU_BOVIN;SAMPLE", "SAMPLE", "crap-ALBU_BOVIN"],
        "Protein.Group": ["cRAP-P02769", "P04264;cRAP-P04264", "cRAP-X1", "cRAP-P02769;X2", "X3", "crap-P02769"],
    }))
    reference = contaminant_reference()
    database = [True, True, False, False, False, False]
    # The twin pair is left to its database category: within a column every member must carry the prefix.
    upstream = [True, False, True, False, False, False]
    both = [True, True, True, False, False, False]
    namespace = {}
    exec(preprocessing_source().code, namespace)
    for accessions in (True, False):
        for categories, expected in (([1, 2, 3, 4, 5], database), ([6], upstream), ([1, 2, 3, 4, 5, 6], both)):
            config = {"enabled": True, "categories": categories, "matchAccessions": accessions}
            assert contaminant_mask(frame, config, reference).tolist() == expected
            assert namespace["contaminant_mask"](frame, config, reference).tolist() == expected
            pd.testing.assert_frame_equal(
                filter_contaminants(frame, enabled=True, categories=categories, match_accessions=accessions),
                frame.loc[[not v for v in expected]])


def test_upstream_category_is_selected_by_default_but_filter_still_needs_enabling():
    assert 6 in ContaminantFilterConfig().categories
    frame = pd.DataFrame({"ProteinNames": ["UNKNOWN"], "Potential contaminant": ["+"]})
    assert not contaminant_mask(frame).any()
    assert contaminant_mask(frame, {"enabled": True}, contaminant_reference()).all()


def test_upstream_category_is_described_separately_in_report():
    from types import SimpleNamespace
    from app.services.report_function import _filter_description, _contaminant_categories_html
    annotation = SimpleNamespace(filter_config=AnnotationFilterConfig(contaminants={"enabled": True, "categories": [6]}),
                                 source_data=pd.DataFrame({"A": [1, 2]}), filtered_data=pd.DataFrame({"A": [2]}))
    assert "independently of cRAP membership" in _filter_description(annotation)
    assert "Cambridge cRAP contaminant filtering" not in _filter_description(annotation)
    assert "Upstream-flagged contaminants" in _contaminant_categories_html(annotation)


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


def test_peptide_reader_preserves_all_rows_without_consulting_protein_filter(monkeypatch, tmp_path):
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
    assert peptide_tools.get_peptide_frame().ProteinNames.tolist() == raw.ProteinNames.tolist()
    assert calls == []
