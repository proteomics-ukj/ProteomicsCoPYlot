"""The plots a user added to the report, and what is written around each one.

An entry is the state of one plot at the moment it was added: its level, the
page's options (`request`), how it was drawn (`render`: Normal or Interactive,
width, height, DPI) and its advanced appearance. Entries are redrawn on the
session's current data when the report is generated.
"""
from __future__ import annotations

import html
import math
from typing import Any, Mapping

import pandas as pd

from app.schemas.stats import AnovaRequest, EnrichmentRequest
from app.schemas.summary import SummaryPlotEntry, SummaryReportContext, SummaryVolcanoEntry
from app.services.comparison_tools import (
    pearson_correlation_html,
    pearson_correlation_png,
    pearson_correlation_table,
    venn_png,
    venn_table,
)
from app.services.method_references import cite
from app.services.plot_styling import migrate_advanced_plot_settings
from app.services.single_protein_tools import (
    single_protein_boxplot_plot,
    single_protein_lineplot_plot,
    unmatched_lineplot_names,
    unmatched_names_notice,
)
from app.services.stats_tools import anova_png, ora_plot_png, run_anova, run_enrichment

# Section key (the text and comment key without its level suffix) and title, in report order.
ENTRY_SECTIONS: tuple[tuple[str, str], ...] = (
    ("Pearson", "Pearson Correlation"),
    ("Venn", "Venn Diagram"),
    ("Volcano", "Volcano Plot"),
    ("Anova", "One-way ANOVA"),
    ("Ora", "GO Enrichment (ORA)"),
    ("Single", "Single Protein"),
)

_CLUSTER_LINKAGE_WORDS = {"average": "average", "complete": "complete", "ward": "Ward", "single": "single"}


def _report():
    from app.services import report_function

    return report_function


def _entries(context: SummaryReportContext, kind: str, *types: str) -> list[SummaryPlotEntry]:
    return [entry for entry in context.entries if entry.kind == kind and entry.type in types]


def _text(request: Mapping[str, Any], key: str, default: str = "") -> str:
    value = request.get(key)
    return default if value is None else str(value).strip()


def _texts(request: Mapping[str, Any], key: str) -> list[str]:
    value = request.get(key)
    if isinstance(value, str):
        value = value.split(",")
    return [str(item).strip() for item in (value or []) if str(item).strip()]


def _flag(request: Mapping[str, Any], key: str, default: bool) -> bool:
    return _report()._to_bool(request.get(key), default)


def _size(entry: SummaryPlotEntry, *, dpi: bool = True) -> dict[str, Any]:
    """The saved figure size, where the entry carries one."""
    bounds = {"widthCm": ("width_cm", 1.0, 60.0), "heightCm": ("height_cm", 1.0, 60.0)}
    if dpi:
        bounds["dpi"] = ("dpi", 72.0, 1200.0)
    size: dict[str, Any] = {}
    for source, (target, low, high) in bounds.items():
        try:
            value = float(entry.render[source])
        except (KeyError, TypeError, ValueError):
            continue
        if not math.isfinite(value):
            continue
        value = min(high, max(low, value))
        size[target] = int(round(value)) if target == "dpi" else value
    return size


def _interactive(entry: SummaryPlotEntry) -> bool:
    return str(entry.render.get("type", "Normal")).strip().lower() == "interactive"


def _paragraph(label: str, text: str) -> str:
    return f"<p><strong>{html.escape(label)}:</strong> {html.escape(text)}</p>"


def _named(count: int, kind: str) -> str:
    singular, plural = _report().entity_labels_for_kind(kind)
    return f"{count:,} {singular if count == 1 else plural}"


def _advanced(entry: SummaryPlotEntry) -> dict[str, Any] | None:
    return None if entry.advanced is None else migrate_advanced_plot_settings(dict(entry.advanced))


def _figure(
    entry: SummaryPlotEntry, warnings: list[str], what: str, png, html_renderer=None, *, args: tuple = (), **options
) -> str:
    report = _report()
    if html_renderer is not None and _interactive(entry):
        height_cm = float(options.get("height_cm", 12.0))
        options.pop("dpi", None)
        return report._safe_plotly(
            renderer=lambda: report.report_plot(html_renderer, *args, **options),
            fallback_title=what,
            warnings=warnings,
            height_px=report._interactive_height_px(height_cm),
        )
    return report._safe_png(
        renderer=lambda: report.report_plot(png, *args, **options),
        fallback_title=what,
        warnings=warnings,
    )


# --------------------------------------------------------------------------- volcano


def volcano_selection(request: Mapping[str, Any]) -> tuple:
    """The settings that decide which comparison a volcano tests."""
    return (
        _text(request, "dataSource", "data"),
        _text(request, "condition1"),
        _text(request, "condition2"),
        _text(request, "testType", "unpaired"),
        float(request.get("pValueThreshold", 0.05)),
        float(request.get("log2fcThreshold", 1.0)),
        bool(request.get("useUncorrected", False)),
    )


def _volcano_entries(
    context: SummaryReportContext, kind: str, warnings: list[str]
) -> tuple[list[SummaryVolcanoEntry], list[dict[str, Any]], list[str]]:
    report = _report()
    entries: list[SummaryVolcanoEntry] = []
    styles: list[dict[str, Any]] = []
    for legacy in context.volcanoEntries:
        if legacy.kind == kind and not legacy.control:
            entries.append(legacy)
            styles.append({})
    for entry in _entries(context, kind, "volcano"):
        fields = {
            key: value for key, value in entry.request.items()
            if key in SummaryVolcanoEntry.model_fields and key not in {"kind", "renderType", "control"}
        }
        try:
            converted = SummaryVolcanoEntry(
                **fields, kind=entry.kind, renderType="Interactive" if _interactive(entry) else "Normal",
            )
        except Exception as exc:
            report._report_failure(warnings, f"{kind.capitalize()} Volcano Plot", exc)
            continue
        entries.append(converted)
        style: dict[str, Any] = _size(entry)
        if entry.advanced is not None:
            style["advanced"] = _advanced(entry)
        styles.append(style)
    anchors = [f"e-{kind}-volcano-{position}" for position in range(1, len(entries) + 1)]
    return entries, styles, anchors


def _volcano_section(context: SummaryReportContext, kind: str, warnings: list[str]) -> tuple[str, dict[tuple, str]]:
    report = _report()
    entries, styles, anchors = _volcano_entries(context, kind, warnings)
    if not entries:
        return "", {}
    fallback = (
        report._global_advanced_settings(context, f"stats:volcano:{kind}") if context.workflowPlotSettings
        else report._global_advanced_settings(context, "stats:volcano")
    )
    content = report._volcano_block_html(
        kind=kind, entries=entries, warnings=warnings, advanced=fallback, styles=styles, anchors=anchors,
    )
    targets: dict[tuple, str] = {}
    for entry, anchor in zip(entries, anchors):
        if f"id='{anchor}'" in content:
            targets.setdefault(volcano_selection(entry.model_dump()), anchor)
    return content, targets


# --------------------------------------------------------------------------- comparison


def _comparison_options(entry: SummaryPlotEntry, names: Mapping[str, str]) -> dict[str, Any]:
    return {target: _text(entry.request, source) for source, target in names.items()}


def _pearson_block(entry: SummaryPlotEntry, index: int, warnings: list[str]) -> str:
    report = _report()
    request, kind = entry.request, entry.kind
    selection = _comparison_options(entry, {
        "mode": "mode", "sample1": "sample1", "sample2": "sample2", "condition1": "condition1",
        "condition2": "condition2", "alias1": "alias1", "alias2": "alias2",
    })
    selection["mode"] = selection["mode"] or "Single"
    what = f"{kind.capitalize()} Pearson Correlation"
    table = pearson_correlation_table(kind=kind, **selection)
    values = dict(zip(table["Metric"], table["Value"]))
    what = f"{what} ({values['Left']} vs {values['Right']})"
    options: dict[str, Any] = dict(
        kind=kind, **selection, color=_text(request, "color"), header=_flag(request, "header", True),
        advanced=_advanced(entry), **_size(entry),
    )
    if "dotSize" in request:
        options["dot_size"] = report._to_float(request.get("dotSize"), 60.0, minimum=0.1, maximum=1000.0)
    figure = _figure(entry, warnings, what, pearson_correlation_png, pearson_correlation_html, **options)
    if not figure:
        return ""
    correlation = values.get("PearsonR", "")
    result = (
        f"Pearson r = {float(correlation):.3f} across {_named(int(values['Points']), kind)} quantified on both sides."
        if correlation != "" else
        f"{_named(int(values['Points']), kind)} quantified on both sides; no correlation could be calculated."
    )
    return (
        f"<h3>Comparison {index}: {html.escape(str(values['Left']))} vs {html.escape(str(values['Right']))}</h3>"
        + report._data_provenance_html(kind=kind, data_source="data", subject="this comparison")
        + _paragraph("Result", result)
        + figure
    )


def _venn_block(entry: SummaryPlotEntry, index: int, warnings: list[str]) -> str:
    report = _report()
    request, kind = entry.request, entry.kind
    selection = _comparison_options(entry, {
        "mode": "mode", "first": "first", "second": "second", "third": "third",
        "alias1": "alias1", "alias2": "alias2", "alias3": "alias3",
    })
    selection["mode"] = selection["mode"] or "Single"
    min_present = report._to_int(request.get("minPresent"), 1, minimum=1)
    labels = [
        selection[alias] or selection[name]
        for alias, name in (("alias1", "first"), ("alias2", "second"), ("alias3", "third")) if selection[name]
    ]
    named = ", ".join(labels[:-1]) + (" and " if len(labels) > 1 else "") + labels[-1] if labels else ""
    what = f"{kind.capitalize()} Venn Diagram ({named})"
    table = venn_table(kind=kind, **selection, min_present=min_present)
    figure = _figure(
        entry, warnings, what, venn_png,
        kind=kind, **selection, min_present=min_present,
        color1=_text(request, "color1"), color2=_text(request, "color2"), color3=_text(request, "color3"),
        header=_flag(request, "header", True), advanced=_advanced(entry), **_size(entry),
    )
    if not figure:
        return ""
    singular, _ = report.entity_labels_for_kind(kind)
    presence = ""
    if selection["mode"].lower() != "single":
        presence = _paragraph(
            "Presence",
            f"A {singular} counts for a condition when it has a value in at least {min_present} of its samples.",
        )
    return (
        f"<h3>Comparison {index}: {html.escape(named)}</h3>"
        + report._data_provenance_html(kind=kind, data_source="data", subject="this comparison")
        + presence
        + figure
        + report._table_to_html("Venn Region Counts", table)
    )


# --------------------------------------------------------------------------- single protein


def _single_block(entry: SummaryPlotEntry, index: int, warnings: list[str]) -> str:
    report = _report()
    request, kind = entry.request, entry.kind
    data_source = report._report_data_source(request.get("dataSource"))
    shared: dict[str, Any] = dict(
        kind=kind, conditions=_texts(request, "conditions"), identifier=_text(request, "identifier", "workflow"),
        data_source=data_source, header=_flag(request, "header", True), legend=_flag(request, "legend", True),
        advanced=_advanced(entry), **_size(entry),
    )
    if entry.type == "single_box":
        protein = _text(request, "protein")
        heading = f"Boxplot {index}: {protein}"
        figure = _figure(
            entry, warnings, f"{kind.capitalize()} Boxplot ({protein})", single_protein_boxplot_plot,
            protein=protein, outliers=_flag(request, "outliers", False), dots=_flag(request, "dots", False),
            sample_n=_flag(request, "sampleN", True),
            random_seed=report._to_int(request.get("randomSeed"), 1, minimum=0, maximum=2_147_483_647),
            **shared,
        )
    else:
        proteins = _texts(request, "proteins")
        shown = ", ".join(proteins[:6]) + (f" and {len(proteins) - 6} more" if len(proteins) > 6 else "")
        heading = f"Line plot {index}: {shown}"
        what = f"{kind.capitalize()} Line Plot ({shown})"
        figure = _figure(
            entry, warnings, what, single_protein_lineplot_plot,
            proteins=proteins, include_id=_flag(request, "includeId", False),
            show_na_low=_flag(request, "showNaLow", False), **shared,
        )
        if figure:
            missing = unmatched_names_notice(unmatched_lineplot_names(
                kind, proteins, identifier=shared["identifier"], data_source=data_source))
            if missing:
                warnings.append(f"{what}: {missing}")
    if not figure:
        return ""
    return (
        f"<h3>{html.escape(heading)}</h3>"
        + report._data_provenance_html(kind=kind, data_source=data_source, subject="this plot")
        + figure
    )


# --------------------------------------------------------------------------- ANOVA


def _anova_request(entry: SummaryPlotEntry) -> tuple[AnovaRequest, int | None]:
    fields = {key: value for key, value in entry.request.items() if key in AnovaRequest.model_fields}
    payload = AnovaRequest.model_validate({**fields, "kind": entry.kind})
    selected = entry.request.get("selectedCluster")
    cluster = None
    if payload.source == "anova_cluster" and str(selected).strip().lower() not in {"", "none", "null"}:
        cluster = int(float(str(selected)))
    return payload, cluster


def anova_cluster_selection(request: Mapping[str, Any]) -> tuple:
    """The settings that decide which rows form which ANOVA cluster."""
    return (
        _text(request, "dataSource", "data"),
        float(request.get("pValueThreshold", 0.05)),
        bool(request.get("useUncorrected", False)),
        int(float(request.get("nClusters", 4))),
        _text(request, "clusterMethod", "average"),
    )


def _anova_statistics(payload: AnovaRequest, result) -> str:
    correction = (
        "uncorrected p-values (no multiple-testing correction)" if payload.useUncorrected
        else f"Benjamini-Hochberg adjusted p-values ({cite('benjamini_hochberg')})"
    )
    threshold_name = "p" if payload.useUncorrected else "adjusted p"
    text = (
        "one-way ANOVA across all conditions (independent groups, equal variances); "
        f"significance at {threshold_name} < {payload.pValueThreshold:g}; {correction}."
    )
    if result.postHocMethod:
        pair = f" for {result.postHocComparison}" if result.postHocComparison else " for every pair of conditions"
        text += f" Tukey HSD post-hoc test{pair} ({cite('tukey_kramer')})."
    return text


def _anova_block(entry: SummaryPlotEntry, index: int, warnings: list[str], anchor: str) -> str:
    report = _report()
    kind = entry.kind
    payload, cluster = _anova_request(entry)
    clustered = payload.source == "anova_cluster"
    if cluster is not None:
        label = f"Cluster {cluster} of {payload.nClusters}"
    else:
        label = "One-way ANOVA, clusters" if clustered else "One-way ANOVA"
    what = f"{kind.capitalize()} ANOVA ({label})"
    figure = _figure(
        entry, warnings, what, anova_png, args=(payload,), selected_cluster=cluster, advanced=_advanced(entry),
    )
    if not figure:
        return ""
    result = run_anova(payload)
    parts = [
        f"<h3 id='{anchor}'>Analysis {index}: {html.escape(label)}</h3>",
        report._data_provenance_html(kind=kind, data_source=payload.dataSource, subject="this analysis", tested=True),
        _paragraph("Statistics", _anova_statistics(payload, result)),
        _paragraph(
            "Result",
            f"{_named(int(result.significantCount), kind)} significant, out of {int(result.totalRows):,} tested."
            + (f" {int(result.postHocSignificantCount):,} significant Tukey HSD comparisons."
               if result.postHocRows else ""),
        ),
    ]
    sizes = pd.DataFrame()
    if clustered:
        linkage = _CLUSTER_LINKAGE_WORDS.get(str(payload.clusterMethod), str(payload.clusterMethod))
        parts.append(_paragraph(
            "Clustering",
            "significant rows, row z-scores, hierarchical clustering on Euclidean distance with "
            f"{linkage} linkage, cut into {payload.nClusters} clusters.",
        ))
        if result.completenessNotice:
            parts.append(_paragraph("Features used", str(result.completenessNotice)))
        if result.clusterRows:
            counts = pd.DataFrame(result.clusterRows)["cluster"].value_counts().sort_index()
            sizes = pd.DataFrame({"Cluster": counts.index.astype(int), "Features": counts.to_numpy()})
    parts.append(figure)
    if cluster is None and not sizes.empty:
        parts.append(report._table_to_html("Cluster Sizes", sizes))
    return "".join(parts)


# --------------------------------------------------------------------------- ORA


_ORA_DIRECTION_WORDS = {"up": "upregulated", "down": "downregulated"}


def _ora_request(entry: SummaryPlotEntry) -> tuple[EnrichmentRequest, str]:
    fields = {key: value for key, value in entry.request.items() if key in EnrichmentRequest.model_fields}
    payload = EnrichmentRequest.model_validate({**fields, "kind": entry.kind})
    direction = _text(entry.request, "direction", "cluster" if payload.source == "anova" else "up")
    if direction not in {"up", "down", "cluster"}:
        raise ValueError("Unknown ORA direction.")
    return payload, direction


def _ora_block(
    entry: SummaryPlotEntry, index: int, warnings: list[str],
    volcano_targets: Mapping[tuple, str], anova_targets: Mapping[tuple, str],
) -> str:
    report = _report()
    kind = entry.kind
    payload, direction = _ora_request(entry)
    from_cluster = direction == "cluster"
    if from_cluster:
        label = f"ANOVA cluster {payload.selectedCluster} of {payload.nClusters}"
    else:
        label = f"{payload.condition1} vs {payload.condition2}, {_ORA_DIRECTION_WORDS[direction]}"
    what = f"{kind.capitalize()} ORA ({label})"
    figure = _figure(
        entry, warnings, what, ora_plot_png, args=(payload,), direction=direction,
        advanced=_advanced(entry), **_size(entry),
    )
    if not figure:
        return ""
    result = run_enrichment(payload)
    genes, terms = {
        "up": (result.upGenes, result.upTerms),
        "down": (result.downGenes, result.downTerms),
        "cluster": (result.clusterGenes, result.clusterTerms),
    }[direction]

    request = payload.model_dump()
    if from_cluster:
        linkage = _CLUSTER_LINKAGE_WORDS.get(str(payload.clusterMethod), str(payload.clusterMethod))
        source = (
            f"cluster {payload.selectedCluster} of {payload.nClusters} from the one-way ANOVA "
            f"({'p' if payload.useUncorrected else 'adjusted p'} < {payload.pValueThreshold:g}, {linkage} linkage)"
        )
        target = anova_targets.get(anova_cluster_selection(request))
    else:
        test = report._TEST_TYPE_REPORT_LABELS.get(str(payload.testType), str(payload.testType))
        controls = (
            f", corrected for {payload.condition1Control} and {payload.condition2Control}"
            if payload.source == "volcano_control" else ""
        )
        source = (
            f"{_ORA_DIRECTION_WORDS[direction]} in the comparison {payload.condition1} vs {payload.condition2}{controls} "
            f"({test}; {'p' if payload.useUncorrected else 'adjusted p'} < {payload.pValueThreshold:g}, "
            f"|log2 fold change| > {payload.log2fcThreshold:g})"
        )
        target = volcano_targets.get(volcano_selection(request)) if payload.source == "volcano" else None
    source_html = f"<p><strong>Source:</strong> {html.escape(source)}; {len(genes):,} genes."
    if target:
        source_html += f" <a href='#{target}'>Go to the underlying plot.</a>"
    source_html += "</p>"

    background = (
        "the genes of this dataset" if str(result.backgroundMode) == "experiment" else "all annotated genes of the database"
    )
    enrichment = (
        f"over-representation of Gene Ontology terms ({cite('gene_ontology')}; {result.organism}) by a one-sided "
        f"hypergeometric test ({cite('hypergeometric_ora')}); background: {background} "
        f"({int(result.backgroundSize):,} genes, {int(result.annotatedBackgroundSize):,} annotated); "
        f"terms with {payload.minTermSize} to {payload.maxTermSize} genes; "
        f"{result.correctionMethod} adjusted p < {payload.adjPValueCutoff:g}; top {payload.topN} terms shown."
    )
    table = pd.DataFrame([
        {
            "Term": term.name,
            "Ontology": term.source,
            "Genes in term": int(term.termSize),
            "Overlap": int(term.intersectionSize),
            "Adjusted p": f"{float(term.adjPValue):.2e}",
        }
        for term in terms
    ])
    return (
        f"<h3>Enrichment {index}: {html.escape(label)}</h3>"
        + source_html
        + report._data_provenance_html(kind=kind, data_source=payload.dataSource, subject="this analysis", tested=True)
        + _paragraph("Enrichment", enrichment)
        + figure
        + (report._table_to_html("Enriched Terms", table) if not table.empty else "")
    )


# --------------------------------------------------------------------------- sections


def _blocks(entries: list[SummaryPlotEntry], warnings: list[str], block, what: str) -> str:
    parts: list[str] = []
    for entry in entries:
        try:
            content = block(entry, len(parts) + 1, warnings)
        except Exception as exc:
            _report()._report_failure(warnings, f"{entry.kind.capitalize()} {what}", exc)
            continue
        if content:
            parts.append(content)
    return "".join(parts)


def optional_sections(
    kind: str, context: SummaryReportContext, warnings: list[str]
) -> list[tuple[str, str, str]]:
    """(section key, title, content) of each section of added plots, in report order."""
    report = _report()
    volcano, volcano_targets = _volcano_section(context, kind, warnings)

    anova_parts: list[str] = []
    anova_targets: dict[tuple, str] = {}
    for position, entry in enumerate(_entries(context, kind, "anova"), start=1):
        anchor = f"e-{kind}-anova-{position}"
        try:
            block = _anova_block(entry, len(anova_parts) + 1, warnings, anchor)
        except Exception as exc:
            report._report_failure(warnings, f"{kind.capitalize()} ANOVA", exc)
            continue
        if not block:
            continue
        anova_parts.append(block)
        if _text(entry.request, "source") == "anova_cluster":
            anova_targets.setdefault(anova_cluster_selection(entry.request), anchor)

    def ora_block(entry: SummaryPlotEntry, index: int, entry_warnings: list[str]) -> str:
        return _ora_block(entry, index, entry_warnings, volcano_targets, anova_targets)

    sections = {
        "Pearson": _blocks(_entries(context, kind, "pearson"), warnings, _pearson_block, "Pearson Correlation"),
        "Venn": _blocks(_entries(context, kind, "venn"), warnings, _venn_block, "Venn Diagram"),
        "Volcano": volcano,
        "Anova": "".join(anova_parts),
        "Ora": _blocks(_entries(context, kind, "ora"), warnings, ora_block, "ORA"),
        "Single": (
            _blocks(_entries(context, kind, "single_box"), warnings, _single_block, "Boxplot")
            + _blocks(_entries(context, kind, "single_line"), warnings, _single_block, "Line Plot")
        ),
    }
    singular, _ = report.entity_labels_for_kind(kind)
    titles = dict(ENTRY_SECTIONS)
    titles["Single"] = f"Single {singular.capitalize()}"
    return [(key, titles[key], sections[key]) for key, _ in ENTRY_SECTIONS]
