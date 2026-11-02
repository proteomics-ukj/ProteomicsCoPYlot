import io
import json
import subprocess
import sys
from zipfile import ZipFile

from app.services import plot_images
from app.services.annotation_store import get_annotation
from app.services.notebook_archive import build_notebook_archive
from app.services.notebook_imputation import applied_imputation_notebook
from app.services.notebook_snapshots import snapshots
from app.services.session_context import get_session_id
from notebook_replay_support import set_original_path


def test_applied_csv_and_three_figures_match_gui(inputs, tmp_path, monkeypatch):
    from unittest.mock import Mock

    diagnostics = Mock(wraps=plot_images.impute_values_with_diagnostics)
    monkeypatch.setattr(plot_images, "impute_values_with_diagnostics", diagnostics)
    token = applied_imputation_notebook("protein", "together", False, {})
    diagnostics.assert_called_once()
    recipe = snapshots.get(get_session_id(), token)
    _, blob = build_notebook_archive(recipe)
    with ZipFile(io.BytesIO(blob)) as archive:
        archive.extractall(tmp_path)
    notebook = json.loads(next(tmp_path.glob("*.ipynb")).read_text())
    cells = "\n".join("".join(cell["source"]) for cell in notebook["cells"] if cell["cell_type"] == "code")
    replay = set_original_path(cells, tmp_path)
    subprocess.run([sys.executable, "-I", "-c", "import sys;sys.path.insert(0,'.')\n" + replay],
                   cwd=tmp_path, check=True, timeout=90)
    annotation = get_annotation("protein")
    assert (tmp_path / "reproduced_imputed.csv").read_text() == annotation.imputed_data.to_csv(index=False).replace("\r\n", "\n")
    options = dict(kind="protein", q_value=0.01, adjust_std=0.3, seed=77, sample_wise=True, advanced={})
    for name, renderer, extra in (("before", plot_images.imputation_before_plot, {}),
                                  ("fit", plot_images.imputation_overall_fit_plot, {"show_imputation_draw": False}),
                                  ("after", plot_images.imputation_after_plot, {"normalization": "together"})):
        assert (tmp_path / "outputs" / name / "reproduced.png").read_bytes() == renderer(**options, **extra)
