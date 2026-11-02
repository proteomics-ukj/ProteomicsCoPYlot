"""Export the applied imputation matrix with its three diagnostic figures."""
from functools import lru_cache

from app.services.notebook_archive import NotebookRecipe, dataframe_bytes
from app.services.notebook_preprocessing import capture_notebook
from app.services.notebook_qc import qc_source
from app.services.notebook_source import extract_notebook_source, function_definition, merge_source, module_source


@lru_cache(maxsize=1)
def imputation_source():
    plots = module_source("plot_images.py")
    science = module_source("functions.py")
    names = ("render_imputation_before", "render_imputation_fit", "render_imputation_after")
    code = [extract_notebook_source(science, ("impute_values_with_diagnostics",)).code,
            qc_source("coverage").code,
            "def _get_plt():\n    return get_plt()",
            "_to_png_bytes = to_png_bytes",
            function_definition(plots, "_imputation_histogram_bin_edges"),
            function_definition(plots, "_imputation_density_histogram_spec")]
    code.extend(function_definition(plots, name) for name in names)
    code.append(function_definition(science, "intensity_unit_label"))
    return extract_notebook_source(merge_source(*code), (*names, "impute_values_with_diagnostics"))


def applied_imputation_notebook(kind, normalization, show_imputation_draw, advanced):
    from app.services.annotation_store import get_annotation
    from app.services.data_tools import _best_imputation_source, _get_sample_columns, metadata_for_kind
    from app.services.dataset_store import get_current_dataset
    from app.services import plot_images
    from app.services.runtime_cache import runtime_cache_version
    from app.services.session_context import get_session_id

    revision = runtime_cache_version()
    annotation, dataset = get_annotation(kind), get_current_dataset(kind)
    info = annotation.imputation if annotation is not None else None
    if dataset is None or info is None or not info.applied or annotation.imputed_data is None:
        raise ValueError("Apply imputation before exporting this notebook.")
    source, frame, _ = _best_imputation_source(kind)
    meta = metadata_for_kind(kind, frame, _get_sample_columns(kind, frame))
    required = [info.seed, info.sampleWise]
    required += [info.downShift, info.width] if info.mode == "numerical" else [info.qValue, info.adjustStd]
    if any(value is None for value in required):
        raise ValueError("The saved imputation settings are incomplete. Apply imputation again to export it.")
    settings = dict(method=info.mode, seed=info.seed, sample_wise=info.sampleWise,
                    q=0.01 if info.qValue is None else info.qValue, adj_std=0.3 if info.adjustStd is None else info.adjustStd,
                    down_shift=1.8 if info.downShift is None else info.downShift, width=0.3 if info.width is None else info.width)
    diagnostics = plot_images.impute_values_with_diagnostics(frame, meta["sample"].tolist(), **settings)
    centered = plot_images.phosprot_values_are_centred(kind)
    references = {
        "before/original.png": plot_images.render_imputation_before(diagnostics, centered, advanced),
        "fit/original.png": plot_images.render_imputation_fit(diagnostics, centered, advanced, info.mode, show_imputation_draw),
        "after/original.png": plot_images.render_imputation_after(diagnostics, centered, advanced, info.mode, normalization),
    }
    recipe = NotebookRecipe(
        plot_name="Applied Imputation", source=imputation_source(), original_filename=dataset.filename,
        provenance={"kind": kind, "plot": "applied-imputation",
                    "dataSource": {"log2": "data", "raw": "source_data"}.get(source, source)},
        resources={"input.parquet": dataframe_bytes(frame), "metadata.parquet": dataframe_bytes(meta),
                   "imputed.parquet": dataframe_bytes(annotation.imputed_data), **references},
        settings={"imputation": settings, "plot": dict(centered=centered, method=info.mode, normalization=normalization,
                  show_imputation_draw=show_imputation_draw, advanced=advanced)},
        calculation="from IPython.display import Image, display\n"
            "diagnostics = functions.impute_values_with_diagnostics(data, metadata['sample'].tolist(), **settings['imputation'])\n"
            "table = diagnostics.imputed_data\n"
            "pd.testing.assert_frame_equal(table, read_resource_table(RESOURCES / 'imputed.parquet'), check_exact=True)\n"
            "table.to_csv('reproduced_imputed.csv', index=False)\ndisplay(table)\n"
            "p = settings['plot']\n"
            "args = (diagnostics, p['centered'], p['advanced'])\n"
            "for name, image in [('before', functions.render_imputation_before(*args)),\n"
            "                    ('fit', functions.render_imputation_fit(*args, p['method'], p['show_imputation_draw'])),\n"
            "                    ('after', functions.render_imputation_after(*args, p['method'], p['normalization']))]:\n"
            "    folder = Path('outputs') / name\n    folder.mkdir(parents=True, exist_ok=True)\n"
            "    (folder / 'reproduced.png').write_bytes(image)\n    display(Image(data=image))",
    )
    if revision != runtime_cache_version():
        raise ValueError("The dataset changed during export. Try again.")
    return capture_notebook(get_session_id(), recipe)
