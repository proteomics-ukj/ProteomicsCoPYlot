import pandas as pd

from app.services import comparison_tools, completeness_tools, phospho_tools, single_protein_tools, table_functions
from app.services.annotation_store import get_annotation
from app.services.data_tools import linear_qc_frame
from app.services.stats_tools import _volcano_dataframe
from proteomicscopylot import notebook_functions as pcp
from proteomicscopylot import notebook_tables
from pathlib import Path


def test_packaged_gui_calculations_are_current():
    from tools.export_notebook_tables import render

    assert Path(notebook_tables.__file__).read_text(encoding="utf-8") == render()


def test_packaged_protein_tables_match_gui(inputs):
    annotation = get_annotation("protein")
    data, meta = annotation.log2_data, annotation.metadata
    for package, gui, options in (
        (pcp.qc_coverage_table, table_functions.qc_coverage_table, {"summary": False}),
        (pcp.qc_coverage_table, table_functions.qc_coverage_table, {"summary": True}),
        (pcp.qc_boxplot_table, table_functions.qc_boxplot_table, {"mode": "Mean"}),
        (pcp.qc_boxplot_table, table_functions.qc_boxplot_table, {"mode": "Single"}),
        (pcp.completeness_sample_summary_table, completeness_tools.completeness_sample_summary_table, {}),
        (pcp.completeness_feature_summary_table, completeness_tools.completeness_feature_summary_table, {}),
    ):
        pd.testing.assert_frame_equal(package(data, meta, **options), gui("protein", **options))
    linear = linear_qc_frame("protein")
    pd.testing.assert_frame_equal(pcp.qc_cv_table(linear, meta), table_functions.qc_cv_table("protein"))
    pd.testing.assert_frame_equal(pcp.qc_feature_cvs(linear, meta), table_functions.qc_feature_cvs(linear, meta))
    for package, gui, options in (
        (pcp.pearson_table, comparison_tools.pearson_correlation_table, {"sample1": "A1", "sample2": "B1"}),
        (pcp.venn_table, comparison_tools.venn_table, {"mode": "Condition", "first": "A", "second": "B"}),
        (pcp.single_protein_lineplot_table, single_protein_tools.single_protein_lineplot_table,
         {"proteins": ["P0", "P1"], "conditions": ["A", "B"]}),
        (pcp.single_protein_boxplot_table, single_protein_tools.single_protein_boxplot_table,
         {"protein": "P0", "conditions": ["A", "B"]}),
    ):
        pd.testing.assert_frame_equal(package(data, meta, "protein", **options), gui("protein", **options))
    for paired in (False, True):
        untested = []
        _volcano_dataframe(kind="protein", condition1="A", condition2="B", identifier="workflow",
                          test_type="paired" if paired else "unpaired", use_uncorrected=False,
                          log2fc_threshold=1, p_value_threshold=.05, untested_out=untested)
        pd.testing.assert_frame_equal(pcp.volcano_untested_table(data, meta, "A", "B", paired=paired), untested[0])


def test_packaged_phosphosite_tables_match_gui(phosphosites):
    data, meta, _ = phospho_tools._phospho_context()
    for package, gui in ((pcp.phospho_counts_table, phospho_tools.phosphosite_plot_table),
                         (pcp.phospho_sty_table, phospho_tools.phospho_sty_table),
                         (pcp.phospho_distribution_table, phospho_tools.phospho_distribution_table)):
        pd.testing.assert_frame_equal(package(data), gui())
    for mode in ("Normal", "Summary"):
        pd.testing.assert_frame_equal(pcp.phospho_coverage_table(data, meta, True, mode=mode),
                                      phospho_tools.phospho_coverage_table(mode=mode))
    annotation = get_annotation("phospho")
    options = dict(kind="phospho", protein="P0", conditions=["A", "B"], filter_m1=False)
    pd.testing.assert_frame_equal(pcp.single_protein_heatmap_table(annotation.log2_data, meta, **options),
                                  single_protein_tools.single_protein_heatmap_table(**options))
