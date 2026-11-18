"""Deterministic, self-contained notebook archives built from trusted plot recipes."""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from importlib.metadata import version
import io
import json
from pathlib import Path, PurePosixPath, PureWindowsPath
import platform
import re
import sys
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

import pandas as pd

from app.app_version import APP_VERSION
from app.services.notebook_source import NotebookSource

_DISTRIBUTIONS = {"PIL": "Pillow", "sklearn": "scikit-learn", "Bio": "biopython",
                  "matplotlib_venn": "matplotlib-venn", "umap": "umap-learn"}


def portable_filename(value: str) -> str:
    return PureWindowsPath(PurePosixPath(value).name).name.replace("\x00", "")


def _json(value) -> bytes:
    return json.dumps(value, ensure_ascii=True, sort_keys=True, indent=2, allow_nan=False).encode("utf-8")


def _cell(kind: str, text: str) -> dict:
    result = {"cell_type": kind, "metadata": {}, "source": text.splitlines(keepends=True)}
    if kind == "code":
        result.update(execution_count=None, outputs=[])
    return result


def dataframe_bytes(frame: pd.DataFrame) -> bytes:
    """Parquet preserves floating-point values, nulls, column types and the index."""
    stream = io.BytesIO()
    portable = frame.copy(deep=False)
    portable.attrs = {}  # Input-file paths or private importer state are not scientific columns.
    portable.to_parquet(stream, engine="pyarrow", index=True)
    return stream.getvalue()


def _fonts_used(value) -> dict[str, str]:
    """The face each font choice in these settings resolved to on this server."""
    from app.services.plot_styling import FONT_ELEMENT_KEYS, FONT_KEY, PLOT_FONTS, installed_family
    keys = {FONT_KEY, *FONT_ELEMENT_KEYS.values()}
    found: dict[str, str] = {}

    def walk(item):
        if isinstance(item, dict):
            for key, child in item.items():
                choice = str(child).strip().lower() if key in keys and child is not None else ''
                if choice in PLOT_FONTS:
                    found[choice] = installed_family(PLOT_FONTS[choice][0])
                else:
                    walk(child)
        elif isinstance(item, (list, tuple)):
            for child in item:
                walk(child)

    walk(value)
    return dict(sorted(found.items()))


# Defined once and used on both sides: the app records a fingerprint with these
# functions and the notebook's first cell checks it with the same text.
NOTEBOOK_HELPERS = """def table_fingerprint(frame, index=False):
    # The content of a table: column names and values, whatever file format or
    # number type they were read with.
    import numpy as np
    digest = hashlib.sha256()
    if index:
        digest.update('\\x1f'.join(map(str, frame.index)).encode('utf-8'))
    for position, column in enumerate(frame.columns):
        values = frame.iloc[:, position]
        digest.update(b'\\x00' + str(column).encode('utf-8') + b'\\x00')
        numbers = pd.to_numeric(values, errors='coerce')
        if not pd.api.types.is_bool_dtype(values) and int(numbers.notna().sum()) == int(values.notna().sum()):
            array = numbers.to_numpy(dtype='float64', copy=True)
            array[np.isnan(array)] = np.nan
            digest.update(array.tobytes())
        else:
            text = values.astype(object).where(values.notna(), '').astype(str)
            digest.update(pd.util.hash_pandas_object(text, index=False).to_numpy().tobytes())
    return digest.hexdigest()


def check_fingerprint(frame, expected, what, index=False):
    # A different table is reported, not refused: the last cell decides whether
    # the outputs are still the ones CoPYlot produced.
    if table_fingerprint(frame, index=index) == expected:
        print(f'Fingerprint matches: {what}.')
        return True
    import warnings
    warnings.warn(f'{what[:1].upper() + what[1:]} is not the table CoPYlot used (its content fingerprint '
                  'differs). The notebook continues; the outputs may differ.', stacklevel=2)
    return False


def original_path(value):
    # An empty path means: use what is shipped with this export.
    text = str(value or '').strip()
    return None if text in {'', '.'} or text.startswith('YOUR_ORIGINAL_') else Path(text)


def select_view(frame, columns, rows=None):
    # The rows and columns of the matrix this output was drawn from.
    return frame.loc[frame.index if rows is None else rows, columns]


def place_imputed_values(frame, values):
    # Imputation is not recomputed: the values CoPYlot drew are placed into the
    # cells that were missing. Rows this output does not use are skipped.
    import numpy as np
    result = frame.copy()
    for column, cells in values.groupby('column', sort=False):
        if column not in result.columns:
            continue
        cells = cells[cells['row'].isin(result.index)]
        series = pd.to_numeric(result[column], errors='coerce').astype('float64')
        series.loc[cells['row'].to_numpy()] = cells['value'].to_numpy(dtype='float64')
        result[column] = series
    return result
"""


def _helpers():
    import hashlib
    namespace = {"hashlib": hashlib, "pd": pd, "Path": Path}
    exec(NOTEBOOK_HELPERS, namespace)
    return namespace


_HELPERS = _helpers()
table_fingerprint = _HELPERS["table_fingerprint"]
place_imputed_values = _HELPERS["place_imputed_values"]
select_view = _HELPERS["select_view"]


# The last cell: every figure or text output is compared with the one CoPYlot
# produced, so "it ran" can no longer be mistaken for "it reproduced".
_VERIFY_CELL = """# Compare every output with the one CoPYlot produced.
results = []
for original in sorted(Path('resources').rglob('original.*')):
    folder = original.parent.relative_to('resources')
    output_dir = Path('outputs') / folder if folder.parts else Path('.')
    reproduced = output_dir / ('reproduced' + original.suffix)
    candidates = sorted(output_dir.glob('reproduced*' + original.suffix))
    if not reproduced.exists() and len(candidates) == 1:
        reproduced = candidates[0]
    if not reproduced.exists():
        results.append(('not written', original))
    else:
        results.append(('identical' if reproduced.read_bytes() == original.read_bytes() else 'DIFFERENT', reproduced))
for status, path in results:
    print(f'{status:>11}  {path}')
if not results:
    print('No figure to compare; the tables were checked exactly above.')
elif all(status == 'identical' for status, _ in results):
    print('Every output is byte-identical to the one CoPYlot produced.')
else:
    from matplotlib import font_manager
    installed = {font.name for font in font_manager.fontManager.ttflist}
    missing = sorted({face for face in manifest.get('fonts', {}).values() if face not in installed})
    hint = (f' The export used the font(s) {", ".join(missing)}, which this machine does not have; '
            'install them or the text is drawn in another face.' if missing else '')
    raise AssertionError('Some outputs differ from the ones CoPYlot produced. The numbers were checked '
        'exactly above, so the difference is in the drawing: most often a different font, or package '
        'versions other than those in requirements.txt.' + hint)
"""


@dataclass(frozen=True)
class NotebookRecipe:
    # Only server-owned adapters may supply executable cells or function names.
    plot_name: str
    source: NotebookSource
    calculation: str
    settings: dict
    provenance: dict
    original_filename: str
    resources: dict[str, bytes]
    uses_imputed_input: bool = False
    imputation_recipe: str = ""
    preparation: str = "data = read_resource_table(RESOURCES / 'input.parquet')"
    input_explanation: str = "The exact input matrix is included to preserve the displayed figure's inputs."
    calculation_sections: tuple[tuple[str, str], ...] = ()
    # A report: builds the combined recipe for the two export options.
    assemble: object = None
    # How the input is rebuilt from the original upload (see `replay_preparation`).
    # None: this output has no such replay and its input matrix is always shipped.
    replay: dict | None = None


def _path_line(variable: str, filename: str, required: bool, placeholder: str, what: str) -> str:
    """One original-file path with the expected file's name beside it."""
    expected = f"{json.dumps(filename, ensure_ascii=True)} " if filename else ""
    if required:
        return f"{variable} = Path('{placeholder}')  # {expected}(required: {what})\n"
    return f"{variable} = ''  # {expected}(optional: while empty, the shipped matrix is used)\n"


def _path_lines(recipe: NotebookRecipe, filename: str) -> str:
    """The original-file path(s), always the first thing a user sets."""
    if recipe.provenance.get('inputReplay') != 'original-file':
        name = f"# Original filename: {json.dumps(filename, ensure_ascii=True)}\n" if filename else ""
        return name + "# Original-file paths, where needed, are configured in the preparation cell below.\n"
    replay = recipe.replay or {}
    required = recipe.provenance.get("requiredOriginalFiles", ["data_path"]) if recipe.replay else ["data_path"]
    lines = _path_line("data_path", filename, "data_path" in required, "YOUR_ORIGINAL_DATA_FILE", "your original uploaded file")
    for variable, name, placeholder, what in replay.get("extraPaths", ()):
        lines += _path_line(variable, portable_filename(name), variable in required, placeholder, what)
    protein = (replay.get("protein") or {}).get("recorded")
    if protein:
        lines += _path_line("protein_data_path", portable_filename(protein["filename"]), "protein_data_path" in required,
                            "YOUR_ORIGINAL_PROTEIN_FILE", "the protein file of the same analysis")
    return lines


# Resources that are matrices of the user's data. Each is shipped only when its
# option is on; everything else in an export is small or is an output.
INPUT_MATRIX = "input.parquet"
BASIS_MATRIX = "input_before_imputation.parquet"
IMPUTED_VALUES = "imputed_values.parquet"


def shipped_matrices(replay: dict, ship_matrices: bool, ship_imputed: bool) -> set[str]:
    """Which data matrices travel in the ZIP for these two options."""
    imputed = replay.get("imputed")
    if not imputed:
        return {INPUT_MATRIX} if ship_matrices else set()
    if imputed == "upstream":
        # A table collapsed from imputed phosphosites: given the drawn values it
        # is deterministic, so either option ships the collapsed matrix.
        return {INPUT_MATRIX} if ship_matrices or ship_imputed else set()
    if ship_imputed:
        return {INPUT_MATRIX}
    return {BASIS_MATRIX} if ship_matrices else set()


def side_resource(side: str, name: str) -> str:
    """A second dataset's files carry its name: protein_input.parquet, protein_processing.json, ..."""
    return f"{side}_{name}" if side else name


def replay_preparation(replay: dict, shipped: set[str], report: dict | None = None, side: str = "") -> str:
    """The code that gets an input matrix of this output.

    With a matrix shipped the notebook runs as it is, and switches to the
    original file as soon as its path is filled in at the top. Without one the
    path is required. Either way the file is checked by content fingerprint.

    In a report (`report` names the section, its path entry and where shared
    matrices are stored) a file is read once and a preprocessing is run once,
    however many sections use it. `side` is a second dataset of the same output
    (the protein table of a protein-corrected site comparison).
    """
    function = replay["function"]
    local = f"{side}_" if side else ""
    data = f"{local}data"
    resource = lambda name: f"RESOURCES / '{side_resource(side, name)}'"  # noqa: E731
    translations = (f"{local}processing['upstream'].get('id_translations', [])" if replay.get("collapse")
                    else f"{local}processing.get('id_translations', [])")
    extra = replay.get("callExtra", "")
    call = f"functions.{function}({local}raw, {local}processing_metadata, {local}processing{extra})"
    what = f"the original {side} file" if side else "the original file"
    check = (f"check_fingerprint(functions.replay_id_translations({local}raw, {translations}), "
             f"{local}recorded['sourceFingerprint'], {what!r})")
    lines = [
        f"{local}processing_metadata = read_resource_table({resource('processing_metadata.parquet')})",
        f"{local}processing = json.loads(({resource('processing.json')}).read_text(encoding='utf-8'))",
        f"{local}rows = (read_resource_table({resource('input_rows.parquet')}).index"
        f" if ({resource('input_rows.parquet')}).exists() else None)",
    ]
    custom = replay.get("custom")
    if custom:
        # An input with its own loading (the long peptide report): read, check, then its recorded steps.
        lines = lines[1:]
        path = "data_path" if report is None else f"data_paths[{json.dumps(report['path'])}]"
        where = ({name: resource(name) for name in (INPUT_MATRIX, BASIS_MATRIX)} if report is None else report["matrices"])
        recorded = "recorded_input" if report is None else f"recorded_input[{json.dumps(report['section'])}]"
        read = custom["read"] if report is None else (
            f"shared(('file', str(original), {custom['readKey']!r}), lambda: {custom['read']})")
        checked = "check_fingerprint(raw, recorded['sourceFingerprint'], 'the original file')"
        if report is not None:
            checked = f"shared(('checked', str(original), recorded['sourceFingerprint']), lambda: {checked})"
        from_file = [f"raw = {read}", checked, *custom["body"], "data = select_view(data, recorded['columns'], rows)"]
    elif report is None:
        path = f"{local}data_path"
        where = {name: resource(name) for name in (INPUT_MATRIX, BASIS_MATRIX)}
        recorded = "recorded_input" + (f"[{side!r}]" if side else "")
        from_file = [f"{local}raw = functions.read_original_table({local}original)", check,
                     f"{data} = select_view({call}, {local}recorded['columns'], {local}rows)"]
    else:
        path, where = f"data_paths[{json.dumps(report['path'])}]", report["matrices"]
        recorded = f"recorded_input[{json.dumps(report['section'])}]" + (f"[{side!r}]" if side else "")
        from_file = [
            f"{local}raw = shared(('file', str({local}original)), lambda: functions.read_original_table({local}original))",
            f"shared(('checked', str({local}original), {local}recorded['sourceFingerprint']), lambda: {check})",
            f"{data} = select_view(shared(('replay', str({local}original), {local}recorded['replayKey']), lambda: {call}), "
            f"{local}recorded['columns'], {local}rows)",
        ]
    lines += [f"{local}recorded = {recorded}", f"{local}original = original_path({path})"]
    imputed = replay.get("imputed") is True
    place = f"{data} = functions.place_imputed_values({data}, read_resource_table({resource('imputed_values.parquet')}))"
    if INPUT_MATRIX in shipped:
        lines += [f"if {local}original is None:",
                  "    # The matrix this output was drawn from, as shipped with the export.",
                  f"    {data} = read_resource_table({where[INPUT_MATRIX]})",
                  *("    " + line for line in (custom or {}).get("shipped", ())),
                  "else:", *("    " + line for line in from_file)]
        if imputed:
            lines.append("    " + place)
    elif BASIS_MATRIX in shipped:
        lines += [f"if {local}original is None:",
                  "    # The matrix before imputation, as shipped with the export.",
                  f"    {data} = read_resource_table({where[BASIS_MATRIX]})",
                  "else:", *("    " + line for line in from_file), place]
    else:
        lines += [f"if {local}original is None:",
                  f"    raise ValueError('Set {path.replace(chr(39), chr(34))} at the top to your original file ' + "
                  f"json.dumps({local}recorded.get('filename', '')) + '. This export ships no data matrix.')",
                  *from_file]
        if imputed:
            lines.append(place)
    label = f"the {side} matrix of this output" if side else "the input matrix of this output"
    lines.append(f"check_fingerprint({data}, {local}recorded['inputFingerprint'], {label!r}, index=True)")
    if replay.get("imputation") and not side:
        before = ("functions.replay_preprocessing(raw, processing_metadata, "
                  "{**processing, 'linear_output': False, 'localization_cutoff': None})")
        if report is not None:
            before = f"shared(('before imputation', str(original), recorded['replayKey']), lambda: {before})"
        lines += [
            "",
            "# How the missing values were imputed in CoPYlot: the same call with the recorded settings.",
            "# It is run and compared here. The notebook continues with the values CoPYlot drew.",
            f"imputation_settings = {json.dumps(replay['imputation'], sort_keys=True)}".replace("true", "True").replace("false", "False").replace("null", "None"),
            "if original is None:",
            "    print('Fill in the original-file path at the top to run the imputation itself. The shipped imputed values are used.')",
            "else:",
            f"    before_imputation = {before}",
            "    imputed_again = functions.impute_values_with_diagnostics(",
            "        before_imputation, processing_metadata['sample'].tolist(), **imputation_settings).imputed_data",
        ]
        if replay.get("linear"):
            lines.append("    print('Imputation was run again with the recorded settings (this output shows linear values, "
                         "so the comparison is left to the final check).')")
        else:
            lines += [
                "    same = table_fingerprint(imputed_again.loc[data.index, recorded['columns']], index=True) == recorded['inputFingerprint']",
                "    print('Imputing again here gives ' + ('exactly the values CoPYlot drew.' if same else",
                "          'other values than CoPYlot drew (another package version or random generator). The shipped values are used.'))",
            ]
    return "\n".join(lines) + "\n"


# The line every protein-corrected adapter reads its protein matrix with; an
# export with a protein replay gets `replay_preparation(..., side='protein')` there.
PROTEIN_MATRIX_LINE = "    protein_data = read_resource_table(RESOURCES / 'protein_input.parquet')\n"


def protein_side_code(code: str, replay: dict, ship_matrices: bool, ship_imputed: bool, report: dict | None = None):
    """`code` with the protein matrix rebuilt like the main input; also which protein matrices ship."""
    side = replay.get("protein")
    if not side or PROTEIN_MATRIX_LINE not in code:
        return code, None
    shipped = shipped_matrices(side, ship_matrices, ship_imputed)
    block = replay_preparation(side, shipped, report, side="protein")
    return code.replace(PROTEIN_MATRIX_LINE, "".join("    " + line + "\n" for line in block.splitlines())), shipped


def export_explanation(replay: dict, shipped: set[str]) -> str:
    what = ("The matrix this output was drawn from is shipped, so the notebook runs as it is. Fill in "
            "the original-file path at the top to rebuild it from your original file instead."
            if shipped else
            "No data matrix is shipped: set the original-file path at the top to your original file.")
    imputed = (" Imputation is not repeated: the values CoPYlot drew are shipped and placed into the missing cells."
               if replay.get("imputed") else "")
    steps = ("" if replay.get("custom") else " and the notebook repeats the recorded identifier translations, "
             "scale handling, normalization, filtering and localization selection")
    return what + " The file is checked by a fingerprint of its content" + steps + "." + imputed + replay.get("notes", "")


def imputed_values_note(text: str, replay: dict | None) -> str:
    """The imputation note of an export that carries the drawn values, not a frozen matrix."""
    if not replay or not replay.get("imputed"):
        return text
    where = ("upstream_imputed_values.parquet holds the phosphosite values CoPYlot drew; they are placed before the collapse."
             if replay["imputed"] == "upstream" else
             "imputed_values.parquet holds the values CoPYlot drew; they are placed into the rebuilt matrix.")
    for old in ("Do not rerun random draws: resources/input.parquet contains the exact imputed values used.",
                "resources/input.parquet contains the exact collapsed phosphoprotein matrix used by this plot."):
        text = text.replace(old, "The outputs never use new random draws: " + where)
    return text


def apply_export_options(recipe: NotebookRecipe, ship_matrices: bool = False, ship_imputed: bool = False) -> NotebookRecipe:
    """The recipe as it is packaged for the user's two export options."""
    from dataclasses import replace

    if recipe.assemble is not None:
        return recipe.assemble(ship_matrices, ship_imputed)
    replay = recipe.replay
    if replay is None:
        return recipe
    shipped = shipped_matrices(replay, ship_matrices, ship_imputed)
    calculation, protein_shipped = protein_side_code(recipe.calculation, replay, ship_matrices, ship_imputed)
    kept = set(shipped) | {side_resource("protein", name) for name in (protein_shipped or ())}
    if replay.get("protein") is None:
        kept.add(side_resource("protein", INPUT_MATRIX))  # A protein matrix without a replay stays frozen.
    matrices = {INPUT_MATRIX, BASIS_MATRIX} | {side_resource("protein", name) for name in (INPUT_MATRIX, BASIS_MATRIX)}
    resources = {name: value for name, value in recipe.resources.items() if name not in matrices or name in kept}
    required = [name for name, needed in (("data_path", not shipped), ("protein_data_path", protein_shipped == set())) if needed]
    return replace(recipe, resources=resources, calculation=calculation,
        imputation_recipe=imputed_values_note(recipe.imputation_recipe, replay),
        preparation=replay_preparation(replay, shipped),
        provenance={**recipe.provenance, "shippedMatrices": sorted(kept & set(resources)),
                    "originalFileRequired": not shipped, "requiredOriginalFiles": required},
        input_explanation=export_explanation(replay, shipped))


def _recorded(replay: dict | None) -> dict:
    """What the notebook checks its inputs against; a second dataset under its name."""
    if not replay:
        return {}
    recorded = dict(replay.get("recorded", {}))
    if replay.get("protein"):
        recorded["protein"] = replay["protein"]["recorded"]
    return recorded


def build_notebook_archive(recipe: NotebookRecipe, *, ship_matrices: bool = False, ship_imputed: bool = False) -> tuple[str, bytes]:
    recipe = apply_export_options(recipe, ship_matrices, ship_imputed)
    if recipe.uses_imputed_input and not recipe.imputation_recipe.strip():
        raise ValueError("Imputed input requires its documented reconstruction recipe.")
    for name in recipe.resources:
        path = PurePosixPath(name)
        if not name or path.is_absolute() or ".." in path.parts or "\\" in name or ":" in name:
            raise ValueError("Resources must have safe relative archive names.")
    if "metadata.parquet" not in recipe.resources or (
            recipe.replay is None and recipe.assemble is None and "input.parquet" not in recipe.resources):
        raise ValueError("The exact input matrix and metadata must be captured.")
    for code in (recipe.preparation, recipe.calculation, *(code for _, code in recipe.calculation_sections)):
        compile(code, "notebook", "exec")
    packages = {"pandas", "pyarrow"}
    packages.update(_DISTRIBUTIONS.get(name, name) for name in recipe.source.imports
                    if name not in sys.stdlib_module_names)
    # Notebook hosting/display are not installed in the production server.
    requirements = []
    for package in sorted(packages - {"jupyterlab", "IPython", "ipython"}):
        requirements.append(f"{package}=={version(package)}")
    requirements.extend(["ipython>=8,<10", "jupyterlab>=4,<5"])
    # The pins are the server's exact versions, which need its Python or newer.
    python_minimum = ".".join(platform.python_version_tuple()[:2])
    python_note = (f"Python {python_minimum} or newer is required (exported with Python "
                   f"{platform.python_version()}); the pinned versions do not install on older Python.")
    filename = portable_filename(recipe.original_filename)
    settings = _json(recipe.settings)
    resources = {f"resources/{name}": value for name, value in recipe.resources.items()}
    resources["resources/settings.json"] = settings
    resources["functions.py"] = recipe.source.code.encode("utf-8")
    manifest = {
        "schemaVersion": 1, "toolVersion": APP_VERSION,
        "pythonVersion": platform.python_version(), "platform": platform.system(),
        "originalFilename": filename, "plot": recipe.plot_name,
        "usesImputedInput": recipe.uses_imputed_input, "provenance": recipe.provenance,
        "requirements": requirements, "pythonRequires": f">={python_minimum}",
        "input": recipe.provenance.get("recordedInputs") or _recorded(recipe.replay),
        # A report keeps each section's settings in its own file.
        "fonts": _fonts_used([recipe.settings, *(json.loads(value) for name, value in recipe.resources.items()
                                                 if name.endswith("settings.json"))]),
        "sha256": {name: sha256(value).hexdigest() for name, value in sorted(resources.items())},
    }
    # Values from filenames, labels and controls live in JSON, never executable cells.
    notebook = {
        "nbformat": 4, "nbformat_minor": 4,
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                     "language_info": {"name": "python", "version": platform.python_version()}},
        "cells": [
            _cell("markdown", f"# {recipe.plot_name}\n\nCoPYlot {APP_VERSION}.\n\n"
                  "Extract the whole ZIP, install requirements.txt and run the cells in order. "
                  "The manifest records the exported inputs and implementation. " + python_note
                  + " The last cell compares every output with the one CoPYlot produced.\n\n" + recipe.input_explanation),
            _cell("code", "import sys\n"
                  f"if sys.version_info < {tuple(int(part) for part in python_minimum.split('.'))!r}:\n"
                  f"    raise RuntimeError({python_note!r})\n"
                  "from pathlib import Path\nimport json\nimport hashlib\nimport pandas as pd\n"
                  "import functions\n"
                  "import pyarrow.parquet as pq\n\n"
                  "def read_resource_table(path):\n"
                  "    # Read synchronously: Arrow dataset workers can race interpreter shutdown.\n"
                  "    # Each resource is one file; no dataset scanner or read-ahead is needed.\n"
                  "    with pq.ParquetFile(path, pre_buffer=False) as parquet:\n"
                  "        return parquet.read(use_threads=False).to_pandas(use_threads=False)\n\n"
                  "RESOURCES = Path('resources')\n"
                  + _path_lines(recipe, filename) +
                  "manifest = json.loads(Path('manifest.json').read_text(encoding='utf-8'))\n"
                  "recorded_input = manifest.get('input', {})\n"
                  "verify_original_export = True  # Set False only after intentionally editing exported files.\n"
                  "if verify_original_export:\n"
                  "    for name, expected in manifest['sha256'].items():\n"
                  "        actual = hashlib.sha256(Path(name).read_bytes()).hexdigest()\n"
                  "        if actual != expected:\n"
                  "            raise ValueError(f'Export file changed: {name}')\n"
                  "settings = json.loads((RESOURCES / 'settings.json').read_text(encoding='utf-8'))\n"
                  "metadata = read_resource_table(RESOURCES / 'metadata.parquet')\n\n" + NOTEBOOK_HELPERS),
        ],
    }
    if recipe.uses_imputed_input:
        notebook["cells"].append(_cell("code", "# Imputation was already performed in CoPYlot. Do not draw values again.\n"
            "# The following documents how the captured input was produced:\n" + "\n".join(
                "# " + line for line in recipe.imputation_recipe.splitlines()) + "\n"))
    notebook["cells"].append(_cell("code", recipe.preparation + "\n"))
    if recipe.calculation_sections:
        for heading, code in recipe.calculation_sections:
            notebook["cells"].extend([_cell("markdown", heading), _cell("code", code + "\n")])
    else:
        notebook["cells"].append(_cell("code", recipe.calculation + "\n"))
    notebook["cells"].append(_cell("code", _VERIFY_CELL))
    stem = re.sub(r"[^A-Za-z0-9_-]+", "_", recipe.plot_name).strip("_")[:100] or "Plot"
    license_root = Path(sys._MEIPASS) / "notebook_sources" if getattr(sys, "frozen", False) else Path(__file__).resolve().parents[3]
    files = {**resources, "LICENSE": (license_root / "LICENSE").read_bytes(),
             "NOTICE": (license_root / "NOTICE").read_bytes(), "manifest.json": _json(manifest),
             "requirements.txt": (f"# {python_note}\n" + "\n".join(requirements) + "\n").encode(),
             f"{stem}.ipynb": _json(notebook)}
    stream = io.BytesIO()
    with ZipFile(stream, "w", compression=ZIP_DEFLATED) as archive:
        for name, content in sorted(files.items()):
            info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, content)
    return stem + ".zip", stream.getvalue()
