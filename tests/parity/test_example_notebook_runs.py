from __future__ import annotations

import contextlib
import io
import json
import sys
import warnings
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
NOTEBOOK = ROOT / "examples" / "example.ipynb"


def _code_cells() -> list[str]:
    payload = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(cell["source"])
        for cell in payload["cells"]
        if cell.get("cell_type") == "code" and "".join(cell["source"]).strip()
    ]


def test_the_example_notebook_exists_and_has_code() -> None:
    assert NOTEBOOK.exists(), NOTEBOOK
    cells = _code_cells()
    assert len(cells) > 30, f"only {len(cells)} code cells found; has the notebook changed shape?"


@pytest.mark.slow
def test_every_cell_of_the_example_notebook_executes() -> None:
    import matplotlib

    matplotlib.use("Agg")

    previous = sys.path[:]
    namespace: dict[str, object] = {"__name__": "__main__"}
    failures: list[str] = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for index, source in enumerate(_code_cells()):
                try:
                    with contextlib.redirect_stdout(io.StringIO()):
                        exec(compile(source, f"<example.ipynb cell {index}>", "exec"), namespace)
                except Exception as exc:
                    failures.append(f"cell {index}: {type(exc).__name__}: {exc}")
    finally:
        sys.path[:] = previous

    assert not failures, (
        "the published example notebook no longer runs against the package:\n  "
        + "\n  ".join(failures[:8])
    )
    _compare_example_with_gui(namespace)


def _compare_example_with_gui(notebook) -> None:
    """Replay the tutorial's documented settings through the GUI endpoints."""
    import pandas as pd
    from fastapi.testclient import TestClient
    from app.main import app
    from proteomicscopylot.notebook_functions import example_data_path

    def compare(actual, expected, key, columns):
        # A repeated display label is not a unique feature ID. Compare every
        # row, including multiplicities, rather than dropping duplicate labels.
        def align(frame):
            return frame[[key, *columns]].sort_values([key, columns[0]]).reset_index(drop=True)
        pd.testing.assert_frame_equal(align(actual), align(expected), check_dtype=False,
                                      check_exact=False, rtol=1e-8, atol=1e-10)

    with TestClient(app) as client:
        def post(path, **kwargs):
            response = client.post(path, **kwargs)
            assert response.status_code == 200, response.text
            return response.json()

        for kind, normalization in (("protein", "median"), ("phospho", "none")):
            for path, name in (("/api/datasets/upload", kind),
                               ("/api/annotations/metadata/upload", kind + "_metadata")):
                source = example_data_path(name)
                with source.open("rb") as handle:
                    post(path, data={"kind": kind}, files={"file": (source.name, handle)})
            post("/api/annotations/generate", json={"kind": kind, "isLog2Transformed": True,
                 "filter": {"minPresent": 3, "mode": "per_group", "normalization": normalization}})
        post("/api/annotations/phosprot/aggregate", json={"mode": "robust_site_median", "source": "non_imputed"})
        retained = client.get("/api/annotations/current/phosprot/filtered/download")
        assert retained.status_code == 200, retained.text
        compare(pd.read_csv(io.StringIO(retained.text)), notebook["phosprot_filtered"], "Phosphoprotein",
                notebook["meta_sites"]["sample"].tolist() + ["site_num"])
        for kind, table, key in (("protein", "volcano", "ProteinNames"),
                                 ("phospho", "site_volcano", "PTM_Collapse_key"),
                                 ("phosprot", "phosprot_volcano", "Phosphoprotein")):
            response = post("/api/stats/volcano/run", json={"kind": kind, "condition1": "WT", "condition2": "GRK",
                            "dataSource": "filtered", "testType": "unpaired"})
            compare(pd.DataFrame(response["rows"]), notebook[table], key,
                    ["log2FC", "pval", "adj_pval", "significance"])


EXAMPLE_R = ROOT / "examples" / "example.R"


from _r_runtime import RSCRIPT


def test_the_r_example_exists() -> None:
    assert EXAMPLE_R.exists(), EXAMPLE_R


@pytest.mark.slow
@pytest.mark.requires_r("plotly")  # the example draws volcano plots
def test_the_r_example_runs(tmp_path) -> None:
    import subprocess

    completed = subprocess.run(
        [RSCRIPT, str(EXAMPLE_R)],
        capture_output=True,
        text=True,
        timeout=900,
        cwd=str(tmp_path),
    )
    if completed.returncode != 0:
        detail = [
            line.strip()
            for line in (completed.stderr or completed.stdout or "").splitlines()
            if line.strip() and "wurde unter R Version" not in line
        ]
        pytest.fail(
            "the published R example no longer runs against the package:\n  "
            + "\n  ".join(detail[-8:])
        )
