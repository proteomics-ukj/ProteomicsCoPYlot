from __future__ import annotations

import pandas as pd


def default_peptide_quantity_column(frame: pd.DataFrame) -> str | None:
    # The importer maps Spectronaut FG.Quantity to the DIA-NN header.
    return "Precursor.Quantity" if "Precursor.Quantity" in frame.columns else None


def resolve_peptide_quantity_column(frame: pd.DataFrame, requested: str | None = None) -> str:
    # Old plot URLs and saved reports may still supply a quantity override.
    column = default_peptide_quantity_column(frame)
    if column is None:
        raise ValueError(
            "Peptide and precursor plots require precursor quantity. Include "
            "Precursor.Quantity (DIA-NN) or FG.Quantity (Spectronaut) in the long report."
        )
    return column
