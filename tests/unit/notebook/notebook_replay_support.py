"""Helpers for replaying exported notebooks in tests."""
import json
import subprocess
import sys

from app.services.annotation_store import get_annotation


def notebook_code(notebook):
    return "\n".join(
        "".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code"
    )


def replay_cells(cells, directory, *, timeout):
    return subprocess.run(
        [sys.executable, "-I", "-c", "import sys; sys.path.insert(0, '.');\n" + cells],
        cwd=directory,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def set_original_path(cells, directory):
    manifest = json.loads((directory / 'manifest.json').read_text())
    provenance = manifest['provenance']
    if provenance.get('inputReplay') != 'original-file':
        return cells
    original = directory / 'test-original.parquet'
    kind = 'phospho' if 'upstream' in provenance['preprocessing'] else provenance['kind']
    if kind in {'peptide', 'precursor'}:
        from app.services.peptide_tools import get_peptide_frame
        frame = get_peptide_frame()
    else:
        frame = get_annotation(kind).source_data
    frame.to_parquet(original, index=True)
    if "Path('YOUR_ORIGINAL_PROTEIN_FILE')" in cells:
        get_annotation('protein').source_data.to_parquet(directory / 'test-original-protein.parquet', index=True)
        cells = cells.replace("Path('YOUR_ORIGINAL_PROTEIN_FILE')", "Path('test-original-protein.parquet')")
    return cells.replace("Path('YOUR_ORIGINAL_DATA_FILE')", "Path('test-original.parquet')")
