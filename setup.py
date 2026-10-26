from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py


class BuildPyWithGoDatabase(_build_py):
    """Copy the bundled GO tables, Cambridge snapshot and their notices into wheels."""

    def run(self) -> None:
        super().run()
        source_dir = Path(__file__).parent / "data" / "db"
        target_dir = Path(self.build_lib) / "proteomicscopylot" / "data"
        self.mkpath(str(target_dir))
        for filename in ("all_go_terms_genes.csv", "all_go_terms_genes_mouse.csv",
                         "crap_ccp.fasta", "crap_ccp_metadata.csv", "CAMBRIDGE_NOTICE.md", "CAMBRIDGE_LICENSE.txt", "CC-BY-4.0.txt"):
            self.copy_file(str(source_dir / filename), str(target_dir / filename))
        self.copy_file(str(source_dir / "PACKAGED_README.md"), str(target_dir / "README.md"))


setup(cmdclass={"build_py": BuildPyWithGoDatabase})
