"""The Replicate & Variance Projection (table and interactive figure) as a notebook export."""
from __future__ import annotations

from functools import lru_cache

from app.services.notebook_archive import NotebookRecipe, dataframe_bytes
from app.services.notebook_preprocessing import capture_notebook, cached_notebook_render, notebook_capture_requested
from app.services.notebook_source import (
    function_definition, merge_source, module_source, portable_module, replace_once,
    extract_notebook_source,
)
from app.services.notebook_volcano import volcano_source

IDENTIFIERS = ["GeneNames", "Gene_group", "ProteinNames", "ProteinAccessions", "PTM_Collapse_key", "UPD_seq",
               "Protein_group", "Phosphoprotein", "site_num", "PTM_Collapse_keys"]


@lru_cache(maxsize=1)
def projection_source():
    original = module_source("stats_tools.py")
    table = function_definition(original, "_simulation_dataframe")
    start, end = table.index("def _simulation_dataframe("), table.index("    if frame_override is not None and metadata_override is not None:")
    table = table[:start] + "def project_table(payload, frame_override, metadata_override, protein_stats=None):\n" + table[end:]
    for before, after in (
        ('''    if frame_override is not None and metadata_override is not None:
        source_used, frame, metadata, warnings = "data", frame_override, metadata_override, []
    else:
        source_used, frame, metadata, warnings = _stats_source(
            payload.kind,
            data_source=payload.dataSource,
        )
''', '''    source_used, frame, metadata, warnings = "data", frame_override, metadata_override, []
'''),
        # The protein table is projected by the caller, from the protein matrix of the export.
        ("            protein_stats = protein_projection_stats(payload)",
         "            raise ValueError('The projected protein comparison is required for protein correction.')"),
        ('    protein_correction = bool(getattr(payload, "proteinCorrection", False)) and frame_override is None',
         '    protein_correction = bool(getattr(payload, "proteinCorrection", False))'),
    ):
        table = replace_once(table, before, after)
    figure = function_definition(original, "simulation_html")
    for before, after in (
        ("def simulation_html(payload: SimulationRequest, advanced: dict | None = None) -> str:",
         "def projection_figure(frame, label_column, payload, advanced=None):"),
        ("    _, label_column, frame, _ = _simulation_dataframe(payload)\n", ""),
        ("    return _plotly_html(fig)", "    return fig"),
    ):
        figure = replace_once(figure, before, after)
    figure = figure.replace("_get_plotly()", "get_plotly_go()")
    base = volcano_source(True)
    statistics = extract_notebook_source(module_source("functions.py"), ("without_variation",)).code
    core = module_source("stats_core.py").replace("bh_adjust", "_bh_adjust").replace("satterthwaite", "_satterthwaite")
    code = merge_source(base.code, statistics, core, table, figure)
    return extract_notebook_source(code, ("project_table", "projection_figure", "plotly_html"))


@lru_cache(maxsize=1)
def _module():
    return portable_module("_copylot_portable_projection", projection_source)


def _protein_stats(module, calculation, frame, metadata):
    from types import SimpleNamespace

    request = SimpleNamespace(**{**calculation, "kind": "protein", "identifier": "workflow", "proteinCorrection": False})
    table = module.project_table(request, frame, metadata)[2]
    return table[["ProteinAccessions", "log2FC", "se", "dof", "pval"]].reset_index(drop=True)


@cached_notebook_render
def render_projection(payload, advanced=None) -> tuple[str, str]:
    from types import SimpleNamespace

    from app.services import stats_tools
    from app.services.annotation_store import get_annotation
    from app.services.dataset_store import get_current_dataset
    from app.services.notebook_provenance import input_history
    from app.services.runtime_cache import runtime_cache_version
    from app.services.session_context import get_session_id

    kind = payload.kind
    annotation, dataset = get_annotation(kind), get_current_dataset(kind)
    if (not notebook_capture_requested() or annotation is None or dataset is None
            or kind not in {"protein", "phospho", "phosprot"}):
        return stats_tools.simulation_html(payload, advanced=advanced), ""
    revision = runtime_cache_version()
    source, data, metadata, warnings = stats_tools._stats_source(kind, data_source=payload.dataSource)
    lineage, imputation, uses_imputed, imputation_text = input_history(kind, annotation, dataset, source=source)
    metadata = metadata[[column for column in ("sample", "condition") if column in metadata]].copy(deep=True)
    columns = list(dict.fromkeys(metadata["sample"].astype(str).tolist() + IDENTIFIERS))
    data = data[[column for column in columns if column in data]].copy(deep=True)
    calculation = payload.model_dump()
    module = _module()
    request = SimpleNamespace(**calculation)
    protein_resources, protein_stats = {}, None
    if payload.proteinCorrection:
        # The protein comparison projected like the sites (same n and variance), from
        # the protein matrix of the same data source. The app's own checks run first.
        stats_tools.protein_projection_stats(payload)
        protein_frame, protein_metadata, _ = stats_tools._protein_frame_for_correction(payload.dataSource, None)
        protein_metadata = protein_metadata[[c for c in ("sample", "condition") if c in protein_metadata]].copy(deep=True)
        protein_columns = list(dict.fromkeys(protein_metadata["sample"].astype(str).tolist() + ["ProteinAccessions", *IDENTIFIERS]))
        protein_frame = protein_frame[[c for c in protein_columns if c in protein_frame]].copy(deep=True)
        protein_stats = _protein_stats(module, calculation, protein_frame, protein_metadata)
        protein_resources = {"protein_input.parquet": dataframe_bytes(protein_frame),
                             "protein_metadata.parquet": dataframe_bytes(protein_metadata)}
    _, label, table, calculation_warnings = module.project_table(request, data, metadata, protein_stats)
    figure = module.projection_figure(table, label, request, advanced)
    content = module.plotly_html(figure)
    recipe = NotebookRecipe(
        plot_name="Replicate and Variance Projection", source=projection_source(),
        settings={"calculation": calculation, "advanced": advanced},
        original_filename=dataset.filename,
        resources={"input.parquet": dataframe_bytes(data), "metadata.parquet": dataframe_bytes(metadata),
                   "table.parquet": dataframe_bytes(table), "original.html": content.encode("utf-8"),
                   "figure.json": figure.to_json().encode("utf-8"), **protein_resources},
        uses_imputed_input=uses_imputed, imputation_recipe=imputation_text,
        provenance={"kind": kind, "plot": "projection", "view": "interactive", "sourceUsed": source,
                    "request": calculation, "labelColumn": label, "warnings": [*warnings, *calculation_warnings],
                    "normalization": annotation.filter_config.normalization,
                    "filter": annotation.filter_config.model_dump(),
                    "upstreamCollapse": lineage.model_dump() if lineage else None,
                    "imputation": imputation.model_dump() if imputation else None},
        calculation="from types import SimpleNamespace\nfrom IPython.display import IFrame, display\n"
            "# A what-if projection: the observed fold changes with the chosen sample size and variance.\n"
            "request = SimpleNamespace(**settings['calculation'])\n"
            "protein_stats = None\n"
            "if request.proteinCorrection:\n"
            "    # Each site's change is corrected for its protein's change, projected the same way.\n"
            "    protein_data = read_resource_table(RESOURCES / 'protein_input.parquet')\n"
            "    protein_metadata = read_resource_table(RESOURCES / 'protein_metadata.parquet')\n"
            "    protein_request = SimpleNamespace(**{**settings['calculation'], 'kind': 'protein',\n"
            "                                         'identifier': 'workflow', 'proteinCorrection': False})\n"
            "    protein_table = functions.project_table(protein_request, protein_data, protein_metadata)[2]\n"
            "    protein_stats = protein_table[['ProteinAccessions', 'log2FC', 'se', 'dof', 'pval']].reset_index(drop=True)\n"
            "_, label_column, table, warnings = functions.project_table(request, data, metadata, protein_stats)\n"
            "pd.testing.assert_frame_equal(table, read_resource_table(RESOURCES / 'table.parquet'), check_exact=True)\n"
            "table.to_csv('reproduced_table.csv', index=False)\ndisplay(table)\n"
            "figure = functions.projection_figure(table, label_column, request, settings['advanced'])\n"
            "html = functions.plotly_html(figure)\n"
            "assert json.loads(figure.to_json()) == json.loads((RESOURCES / 'figure.json').read_text())\n"
            "Path('reproduced.html').write_text(html, encoding='utf-8', newline='')\n"
            "display(IFrame('reproduced.html', width='100%', height=760))\n",
    )
    try:
        token = capture_notebook(get_session_id(), recipe) if runtime_cache_version() == revision else ""
    except ValueError:
        token = ""
    return content, token
