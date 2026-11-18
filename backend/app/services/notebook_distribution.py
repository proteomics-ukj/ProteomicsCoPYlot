"""The QQ plot of the Distribution page as a notebook export."""
from __future__ import annotations

from functools import lru_cache

from app.schemas.annotation import AnnotationKind
from app.services.notebook_archive import NotebookRecipe, dataframe_bytes
from app.services.notebook_preprocessing import capture_notebook, cached_notebook_render, notebook_capture_requested
from app.services.notebook_source import (
    bound_settings, extract_notebook_source, function_definition, merge_source, module_source,
    portable_module, portable_output_format, portable_styling, replace_once,
)


@lru_cache(maxsize=1)
def qqnorm_source():
    rendering = module_source("plot_images.py")
    function = function_definition(rendering, "distribution_qqnorm_plot")
    for before, after in (
        ("def distribution_qqnorm_plot(", "def plot_figure("),
        ("    kind: AnnotationKind,", "    data, metadata,"),
        ("plt = _get_plt()", "plt = get_plt()"),
        ("frame = log2_qc_frame(kind)", "frame = data"),
        ("sample_columns = _get_sample_columns(kind, frame)", "sample_columns = metadata['sample'].tolist()"),
        ('raise ValueError(f"No sample columns found for {kind} dataset.")', 'raise ValueError("No annotated sample columns found.")'),
        ("return _to_png_bytes(fig, plt, dpi=300)", "return fig"),
    ):
        function = replace_once(function, before, after)
    helpers = function_definition(rendering, "_make_fig").replace("_cm_to_inch(", "cm_to_inch(")
    science = module_source("functions.py")
    quantiles = "QQ_MAX_POINTS = 100_000\n\n" + extract_notebook_source(science, ("qqnorm_plot_data",)).code
    styles = extract_notebook_source(portable_styling(include_matplotlib=True),
        ("marker_style", "resolve_reference_line", "apply_advanced_axes", "get_plt")).code
    colors = extract_notebook_source(module_source("_core.py"), ("resolve_plot_colors",)).code
    serialization = extract_notebook_source(module_source("plot_export.py"), ("to_png_bytes", "cm_to_inch")).code
    combined = merge_source("import math\nimport numpy as np\nimport pandas as pd\nfrom typing import Any",
                            colors, quantiles, styles, serialization, helpers, function)
    return extract_notebook_source(combined, ("plot_figure", "to_png_bytes", "get_plt"))


@lru_cache(maxsize=1)
def _module():
    return portable_module("_copylot_portable_qqnorm", qqnorm_source)


@cached_notebook_render
def render_qqnorm(kind: AnnotationKind, **options) -> tuple[bytes, str]:
    from app.services import plot_images
    from app.services.annotation_store import get_annotation
    from app.services.data_tools import _get_sample_columns, log2_qc_frame, metadata_for_kind
    from app.services.dataset_store import get_current_dataset
    from app.services.notebook_provenance import input_history
    from app.services.plot_export import plot_output_format
    from app.services.runtime_cache import runtime_cache_version
    from app.services.session_context import get_session_id

    if not notebook_capture_requested() or kind not in {"protein", "phospho", "phosprot"}:
        return plot_images.distribution_qqnorm_plot(kind, **options), ""
    revision = runtime_cache_version()
    annotation, dataset = get_annotation(kind), get_current_dataset(kind)
    if annotation is None or dataset is None:
        return plot_images.distribution_qqnorm_plot(kind, **options), ""
    lineage, imputation, uses_imputed_input, imputation_text = input_history(kind, annotation, dataset)
    data = log2_qc_frame(kind).copy(deep=True)
    metadata = metadata_for_kind(kind, data, _get_sample_columns(kind, data)).copy(deep=True)
    settings = bound_settings(plot_images.distribution_qqnorm_plot, kind, **options)
    settings.pop("kind")
    module = _module()
    figure = module.plot_figure(data, metadata, **settings)
    fmt = plot_output_format.get()
    with portable_output_format(module, fmt):
        content = module.to_png_bytes(figure, module.get_plt(), dpi=300)
    recipe = NotebookRecipe(
        plot_name="QQ Plot", source=qqnorm_source(), settings=settings, original_filename=dataset.filename,
        provenance={"kind": kind, "plot": "qqnorm", "inputScale": "log2", "dataSource": "data",
                    "upstreamCollapse": lineage.model_dump() if lineage else None,
                    "imputation": imputation.model_dump() if imputation else None,
                    "isLog2Transformed": annotation.is_log2_transformed, "outputs": ["plot"],
                    "renderer": "Exact UI QQ plot renderer with portable inputs"},
        resources={"input.parquet": dataframe_bytes(data), "metadata.parquet": dataframe_bytes(metadata),
                   f"original.{fmt}": bytes(content)},
        uses_imputed_input=uses_imputed_input, imputation_recipe=imputation_text,
        calculation="from IPython.display import Image, display\n"
                    "figure = functions.plot_figure(data, metadata, **settings)\n"
                    "image = functions.to_png_bytes(figure, functions.get_plt(), dpi=300)\n"
                    "Path('reproduced.png').write_bytes(image)\ndisplay(Image(data=image))\n",
    )
    try:
        token = capture_notebook(get_session_id(), recipe) if runtime_cache_version() == revision else ""
    except ValueError:
        token = ""
    return content, token
