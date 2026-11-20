"""The Kinase Library volcano, drawn from two uploaded result files, as a notebook export.

Its input is not the dataset: it is the pair of result tables the user got from
the Kinase Library site. They are small and travel in the ZIP as uploaded.
"""
from __future__ import annotations

from functools import lru_cache

import pandas as pd

from app.services.notebook_archive import NotebookRecipe, dataframe_bytes
from app.services.notebook_snapshots import snapshots
from app.services.notebook_source import (
    extract_notebook_source, function_definition, merge_source, module_source, portable_module, portable_styling,
    replace_once,
)


@lru_cache(maxsize=1)
def kinase_volcano_source():
    original = module_source("phospho_tools.py")
    figure = function_definition(original, "ksea_uploaded_volcano_html")
    start, end = figure.index("    if evidence_out is not None:"), figure.index("    color_map = {")
    figure = figure[:start] + figure[end:]
    for before, after in (
        ("def ksea_uploaded_volcano_html(", "def plot_figure("),
        ("    evidence_out: dict | None = None,\n", ""),
        ("go = _get_plotly_go()", "go = get_plotly_go()"),
        ("    return _plotly_html(fig)", "    return fig"),
    ):
        figure = replace_once(figure, before, after)
    constants = ('KINASE_ENRICHED_UP = "Enriched in up-regulated sites"\n'
                 'KINASE_ENRICHED_DOWN = "Enriched in down-regulated sites"\n')
    styles = extract_notebook_source(portable_styling(include_matplotlib=False),
        ("plotly_marker", "resolve_reference_line", "advanced_plotly_layout")).code
    science = extract_notebook_source(module_source("functions.py"), ("comparison_contrast",)).code
    common = module_source("plotting_common.py")
    html = function_definition(common, "plotly_html")
    html = replace_once(html, "full_html=False,", "full_html=False, div_id='copylot-plot',\n        post_script=\"window.dispatchEvent(new Event('copylot-plot-ready'));\",")
    code = merge_source("import re\nimport numpy as np\nimport pandas as pd\nfrom typing import Any", constants, styles, science,
        function_definition(original, "_normalize_ksea_results"), figure, html,
        function_definition(common, "get_plotly"), function_definition(common, "get_plotly_go"),
        extract_notebook_source(module_source("plot_webgl.py"), ("force_svg_rendering",)).code,
        extract_notebook_source(module_source("html_security.py"), ("GENERATED_HTML_META",)).code)
    return extract_notebook_source(code, ("plot_figure", "plotly_html"))


@lru_cache(maxsize=1)
def _module():
    return portable_module("_copylot_portable_kinase_volcano", kinase_volcano_source)


def render_kinase_volcano(*, st_results: pd.DataFrame, tnc_results: pd.DataFrame, st_name: str = "", tnc_name: str = "",
                          **options) -> tuple[str, str]:
    """The figure as HTML and the token of its export ('' when exports are off here)."""
    from app.services.notebook_preprocessing import notebook_capture_requested
    from app.services.session_context import get_session_id

    module = _module()
    figure = module.plot_figure(st_results=st_results, tnc_results=tnc_results, **options)
    content = module.plotly_html(figure)
    if not notebook_capture_requested():
        return content, ""
    empty = dataframe_bytes(pd.DataFrame())
    recipe = NotebookRecipe(
        plot_name="Kinase Library Volcano", source=kinase_volcano_source(), settings=dict(options),
        original_filename="",
        resources={"input.parquet": empty, "metadata.parquet": empty,
                   "st_results.parquet": dataframe_bytes(st_results), "tnc_results.parquet": dataframe_bytes(tnc_results),
                   "original.html": content.encode("utf-8"), "figure.json": figure.to_json().encode("utf-8")},
        provenance={"kind": "phospho", "plot": "kinase-library-volcano", "view": "interactive",
                    "inputReplay": "uploaded-results", "uploadedFiles": {"st": st_name, "tnc": tnc_name},
                    "renderer": "Exact UI Kinase Library volcano renderer"},
        input_explanation="This figure is drawn from the two Kinase Library result files you uploaded "
            "(serine/threonine and tyrosine). They are the input and travel in the ZIP as uploaded "
            "(st_results.parquet, tnc_results.parquet); your dataset is not read.",
        calculation="from IPython.display import IFrame, display\n"
            "st_results = read_resource_table(RESOURCES / 'st_results.parquet')\n"
            "tnc_results = read_resource_table(RESOURCES / 'tnc_results.parquet')\n"
            "figure = functions.plot_figure(st_results=st_results, tnc_results=tnc_results, **settings)\n"
            "html = functions.plotly_html(figure)\n"
            "assert json.loads(figure.to_json()) == json.loads((RESOURCES / 'figure.json').read_text())\n"
            "Path('reproduced.html').write_text(html, encoding='utf-8', newline='')\n"
            "display(IFrame('reproduced.html', width='100%', height=760))\n",
    )
    try:
        return content, snapshots.capture(get_session_id(), recipe)
    except ValueError:
        return content, ""
