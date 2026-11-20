from __future__ import annotations

import numpy as np
import pandas as pd

from app.schemas.annotation import AnnotationKind
from app.services.data_tools import _get_sample_columns, metadata_for_kind as _metadata_for_kind, linear_qc_frame, log2_qc_frame
from app.services.functions import (
    sample_display_labels,
    condition_cvs,
    extract_sample_id as _extract_id_or_number,
    raw_zero_to_nan,
)
from app.services.plot_images import _grouped_meta_order
from app.services.runtime_cache import apply_cached_wrappers


def prepare_coverage_df(data: pd.DataFrame, meta: pd.DataFrame, name: bool = False) -> pd.DataFrame:
    data = data.copy()
    meta = _grouped_meta_order(meta.copy())
    meta["sample"] = meta["sample"].astype(str)
    meta["new_sample"] = meta.groupby("condition").cumcount() + 1
    meta['new_sample'] = sample_display_labels(meta, index_column="new_sample")

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data_filtered = data.rename(columns=rename_dict)
    annotated_columns = [c for c in meta["new_sample"].tolist() if c in data_filtered.columns]
    if not annotated_columns:
        return pd.DataFrame(columns=["Sample", "Number"])

    data_filtered = data_filtered[annotated_columns].notna().astype(int)
    data_long = data_filtered.reset_index(drop=True).melt(var_name="Sample", value_name="Number")

    if name:
        data_annotated = data_long.merge(
            meta[["sample", "new_sample"]],
            left_on="Sample",
            right_on="new_sample",
            how="left",
        )
        summary = (
            data_annotated.groupby(["Sample", "sample"], as_index=False)["Number"]
            .sum()
            .rename(columns={"sample": "Original_Name"})
        )
        summary = summary[["Sample", "Number", "Original_Name"]]
    else:
        summary = data_long.groupby("Sample", as_index=False)["Number"].sum()

    summary["Sample"] = pd.Categorical(summary["Sample"], categories=annotated_columns, ordered=True)
    summary = summary.sort_values(["Sample"])

    # The union of the per-sample counts above, not the row count of the matrix.
    # A feature can be identified -- real sequence, localization, protein -- with
    # no quantity surviving in any run, and those rows are in the matrix but in
    # none of the bars. Counting them here made the total larger than the union
    # of the rows beneath it, and disagree with the missing-value plot, which
    # already excludes them.
    quantified_anywhere = data_filtered.sum(axis=1) > 0
    total_value = int(quantified_anywhere.sum())

    all_samples_row: dict[str, object] = {"Sample": "All samples", "Number": total_value}
    if name:
        all_samples_row["Original_Name"] = ""
    extra_rows = [all_samples_row]

    # Never silently dropped: an identified-but-unquantified feature is a real
    # finding (two of the HEK293 ones are Class I), so the count is stated rather
    # than folded into the total or left out of the table altogether.
    identified_only = int(len(quantified_anywhere) - total_value)
    if identified_only:
        unquantified_row: dict[str, object] = {
            "Sample": "Identified, not quantified",
            "Number": identified_only,
        }
        if name:
            unquantified_row["Original_Name"] = ""
        extra_rows.append(unquantified_row)

    summary = pd.concat([summary, pd.DataFrame(extra_rows)], ignore_index=True)
    return summary


def prepare_coverage_summary_df(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    data = data.copy()
    meta = _grouped_meta_order(meta.copy())
    meta["sample"] = meta["sample"].astype(str)
    meta["condition"] = meta["condition"].astype(str)
    condition_order = meta["condition"].dropna().unique().tolist()
    meta["id"] = meta["sample"].apply(_extract_id_or_number)
    meta["new_sample"] = meta.groupby("condition").cumcount() + 1
    meta['new_sample'] = sample_display_labels(meta, index_column="new_sample", include_id=True)

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data_filtered = data.rename(columns=rename_dict)
    annotated_cols = [c for c in meta["new_sample"].tolist() if c in data_filtered.columns]
    if not annotated_cols:
        return pd.DataFrame(columns=["condition", "mean", "median", "sd", "min", "max"])

    data_binary = data_filtered[annotated_cols].notna().astype(int)
    # Count by column so samples with no remaining feature rows still have zero coverage.
    sample_summary = data_binary.sum().rename_axis("Sample").reset_index(name="Value")
    sample_summary = sample_summary.merge(
        meta[["new_sample", "condition"]], left_on="Sample", right_on="new_sample"
    )
    condition_summary = (
        sample_summary.groupby("condition", sort=False)
        .agg(mean=("Value", "mean"), median=("Value", "median"), sd=("Value", "std"),
             min=("Value", "min"), max=("Value", "max"))
        .reset_index()
    )

    if condition_order:
        condition_summary["condition"] = pd.Categorical(
            condition_summary["condition"],
            categories=condition_order,
            ordered=True,
        )
        condition_summary = condition_summary.sort_values("condition").reset_index(drop=True)
        condition_summary["condition"] = condition_summary["condition"].astype(str)
    # The same figures over every sample, whatever its condition.
    counts = sample_summary["Value"]
    overall = pd.DataFrame([{
        "condition": "All samples", "mean": counts.mean(), "median": counts.median(), "sd": counts.std(),
        "min": counts.min(), "max": counts.max(),
    }])
    condition_summary = pd.concat([condition_summary, overall], ignore_index=True)
    rounded = ["mean", "median", "sd"]
    condition_summary[rounded] = condition_summary[rounded].astype(float).round(2)
    condition_summary[["min", "max"]] = condition_summary[["min", "max"]].astype(int)
    return condition_summary


def _boxplot_summary(values):
    from matplotlib.cbook import boxplot_stats

    stats = boxplot_stats(values, whis=1.5)[0]
    return {
        name: round(float(stats[key]), 2)
        for name, key in (
            ("mean", "mean"), ("median", "med"), ("Q1", "q1"), ("Q3", "q3"),
            ("lower_whisker", "whislo"), ("upper_whisker", "whishi"),
        )
    }


def prepare_boxplot_df(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    conditions_order = meta["condition"].unique()
    summary_list: list[dict[str, object]] = []
    for condition in conditions_order:
        columns = meta.loc[meta["condition"] == condition, "sample"].tolist()
        values = data[columns].mean(axis=1, skipna=True).values
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue

        summary = _boxplot_summary(values)

        summary_list.append(
            {
                "Condition": str(condition),
                **summary,
            }
        )

    summary_df = pd.DataFrame(summary_list)
    if summary_df.empty:
        return summary_df
    return summary_df.reset_index(drop=True)


def prepare_boxplot_single_df(data: pd.DataFrame, meta: pd.DataFrame, name: bool = False) -> pd.DataFrame:
    meta = meta.copy()
    meta["id"] = meta["sample"].apply(_extract_id_or_number)
    meta["new_sample"] = meta.groupby("condition").cumcount() + 1
    meta['new_sample'] = sample_display_labels(meta, index_column="new_sample")

    rename_dict = dict(zip(meta["sample"], meta["new_sample"]))
    data_renamed = data.rename(columns=rename_dict)
    samples = [s for s in meta["new_sample"].tolist() if s in data_renamed.columns]

    summary_list: list[dict[str, object]] = []
    for sample in samples:
        values = data_renamed[sample].values
        values = values[np.isfinite(values)]
        if len(values) == 0:
            continue

        summary = _boxplot_summary(values)

        row_dict: dict[str, object] = {
            "Sample": sample,
            **summary,
        }
        if name:
            original_name = meta.loc[meta["new_sample"] == sample, "sample"].values[0]
            row_dict["Original_Name"] = str(original_name)
        summary_list.append(row_dict)

    summary_df = pd.DataFrame(summary_list)
    if summary_df.empty:
        return summary_df
    return summary_df.reset_index(drop=True)


def prepare_cov_df(data: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    summary_list: list[dict[str, object]] = [
        {"Condition": condition, **_boxplot_summary(cv)} for condition, cv in condition_cvs(data, meta)
    ]

    summary_df = pd.DataFrame(summary_list)
    if summary_df.empty:
        return summary_df
    return summary_df.reset_index(drop=True)


def qc_coverage_table(kind: AnnotationKind, summary: bool = False) -> pd.DataFrame:
    frame = log2_qc_frame(kind)
    meta = _metadata_for_kind(kind, frame, _get_sample_columns(kind, frame))
    if summary:
        return prepare_coverage_summary_df(frame, meta)
    return prepare_coverage_df(frame, meta)


def qc_boxplot_table(kind: AnnotationKind, mode: str = "Mean") -> pd.DataFrame:
    frame = log2_qc_frame(kind)
    sample_columns = _get_sample_columns(kind, frame)
    meta = _metadata_for_kind(kind, frame, sample_columns)
    numeric = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    if str(mode).lower() == "single":
        return prepare_boxplot_single_df(numeric, meta)
    return prepare_boxplot_df(numeric, meta)


def qc_feature_cvs(frame: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Per feature: mean log2 quantity and CV (%) of every condition, then of all samples.

    Everything is calculated here from the linear sample values; of the uploaded
    table only the identifiers are kept (names, genes, accessions, as far as the
    table has them), so the columns are the same whichever software wrote it.
    """
    samples = meta["sample"].tolist()
    frame = frame.reset_index(drop=True)
    values = raw_zero_to_nan(frame[samples], samples)
    cvs = dict(condition_cvs(values, meta))
    with np.errstate(divide="ignore", invalid="ignore"):
        log2_values = np.log2(values.where(values > 0))
    identifiers = ("ProteinNames", "GeneNames", "ProteinAccessions",
                   "PTM_Collapse_key", "Phosphoprotein", "Protein_group", "Gene_group")
    result = frame[[column for column in identifiers if column in frame.columns]].copy()
    if result.shape[1] == 0:
        result["Feature row"] = frame.index + 1

    def add(label: str, column) -> None:
        while label in result.columns:
            label += " (calculated)"
        result[label] = column

    conditions = meta["condition"].astype(str)
    for condition in conditions.unique():
        columns = meta.loc[conditions == condition, "sample"].tolist()
        add(f"{condition} mean log2", log2_values[columns].mean(axis=1, skipna=True))
        add(f"{condition} CV (%)", cvs.get(condition, np.nan))
    add("All samples mean log2", log2_values[samples].mean(axis=1, skipna=True))
    overall = values[samples].std(axis=1, skipna=True) / values[samples].mean(axis=1, skipna=True) * 100.0
    add("All samples CV (%)", overall.where(np.isfinite(overall)))
    return result


def qc_cv_table(kind: AnnotationKind) -> pd.DataFrame:
    frame = linear_qc_frame(kind)
    sample_columns = _get_sample_columns(kind, frame)
    frame = raw_zero_to_nan(frame, sample_columns)
    meta = _metadata_for_kind(kind, frame, sample_columns)
    numeric = frame[sample_columns].apply(pd.to_numeric, errors="coerce")
    return prepare_cov_df(numeric, meta)


apply_cached_wrappers(
    globals(),
    [
        "qc_coverage_table",
        "qc_boxplot_table",
        "qc_cv_table",
    ],
)
