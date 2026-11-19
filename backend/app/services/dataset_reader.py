from __future__ import annotations

import os
import re
from typing import BinaryIO, Callable

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from app.services.upload_limits import (
    ensure_columns_within_budget,
    ensure_frame_bytes_within_budget,
    ensure_table_within_cell_budget,
)


ALLOWED_EXTENSIONS = {"csv", "tsv", "txt", "xlsx", "parquet"}


def get_extension(filename: str) -> str:
    ext = os.path.splitext(filename)[1].lower().lstrip(".")
    if ext not in ALLOWED_EXTENSIONS:
        raise ValueError(f"Unsupported file type: {ext}")
    return ext


def make_columns_unique(columns: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    result: list[str] = []
    occupied = {str(col) for col in columns}

    for col in columns:
        base = str(col)
        if base not in seen:
            seen[base] = 0
            result.append(base)
        else:
            seen[base] += 1
            candidate = f"{base}.{seen[base]}"
            while candidate in occupied:
                seen[base] += 1
                candidate = f"{base}.{seen[base]}"
            occupied.add(candidate)
            result.append(candidate)
    return result


def _normalize_colname(name: object) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).strip().lower())


def rename_cols(df: pd.DataFrame) -> pd.DataFrame:
    # Trim whitespace first to make downstream aliasing predictable.
    renamed = df.copy()
    renamed.columns = [str(col).strip() for col in renamed.columns]

    alias_targets = {
        "run": "File.Name",
        "rfilename": "File.Name",
        # The accession column, kept apart from ProteinNames: Spectronaut's
        # PG.ProteinNames holds entry names (TMA7_HUMAN) and the accessions sit
        # in PG.ProteinAccessions; DIA-NN's Protein.Group, FragPipe's "Protein
        # ID" and MaxQuant's "Majority protein IDs" are accessions. Phosphosite
        # rows name their protein by accession (Protein_group), so this is the
        # column a site is matched to its protein by.
        "pgproteinaccessions": "ProteinAccessions",
        "proteingroup": "ProteinAccessions",
        "proteinid": "ProteinAccessions",
        "majorityproteinids": "ProteinAccessions",
        "proteinaccessions": "ProteinAccessions",
        "pgproteinnames": "ProteinNames",
        "pgproteingroups": "ProteinNames",
        "proteinids": "ProteinNames",
        "proteinnames": "ProteinNames",
        # FragPipe combined_protein.tsv: "Protein" (the full FASTA header) and
        # "Gene". "Protein ID" (the bare accession) stays as it is.
        "protein": "ProteinNames",
        "genes": "GeneNames",
        "genenames": "GeneNames",
        "pggenes": "GeneNames",
        "gene": "GeneNames",
        "pepstrippedsequence": "Stripped.Sequence",
        "egmodifiedsequence": "Modified.Sequence",
        "fgquantity": "Precursor.Quantity",
        "egrtpredicted": "Predicted.RT",
        "egapexrt": "RT",
    }

    current_columns = [str(col) for col in renamed.columns]
    occupied = set(current_columns)
    rename_map: dict[str, str] = {}
    for column in current_columns:
        # `Protein_group` is CoPYlot's own phosphosite column (the site's
        # protein accession, from the collapse); it normalises like DIA-NN's
        # "Protein.Group" but must keep its name.
        if column == "Protein_group":
            continue
        target = alias_targets.get(_normalize_colname(column))
        if not target or target == column:
            continue
        # Avoid creating unnecessary duplicate canonical columns.
        if target in occupied:
            continue
        rename_map[column] = target
        occupied.add(target)
    renamed = renamed.rename(columns=rename_map)
    if ("Precursor.Id" not in renamed and "EG.PrecursorId" not in renamed
            and {"FG.LabeledSequence", "FG.Charge"}.issubset(renamed.columns)):
        sequence = renamed["FG.LabeledSequence"].astype("string").str.strip()
        charge = pd.to_numeric(renamed["FG.Charge"], errors="coerce")
        charge = charge.where(charge.gt(0) & charge.mod(1).eq(0)).astype("Int64")
        renamed["Precursor.Id"] = sequence.where(sequence.ne("")) + "/" + charge.astype("string")

    # DIA-NN's Protein.Ids is both the label column (it keeps feeding
    # ProteinNames, as it always has) and an accession list; when the file has
    # no separate accession column, it serves as ProteinAccessions too.
    if "ProteinAccessions" not in renamed.columns:
        for column in current_columns:
            if _normalize_colname(column) == "proteinids" and column in rename_map:
                renamed["ProteinAccessions"] = renamed[rename_map[column]]
                break
    return renamed


# Reports of other tools ------------------------------------------------------
#
# `adapt_report` brings six exports into the shapes the app reads: a matrix
# (one row per protein or site, one column per run) for the table datasets and
# a long report (one row per precursor and run) for the peptide dataset. It
# only reshapes and renames; values are not transformed. Tables it does not
# recognise pass through unchanged. Decided 2026-11-13 on PRIDE examples
# (PXD071383, PXD071048, PXD085265, PXD077951), see D34; the FragPipe peptide
# tables were added 2026-11-18 on PXD076585, see D37.

# "[1] GFP DMSO1.PG.Quantity", "255.raw.PTM.Quantity", "[3] run.raw.EG.TotalQuantity (Settings)"
_SPECTRONAUT_RUN_COLUMN = re.compile(r"^(?P<run>.+)\.(?P<level>PG|PEP|EG|FG|PTM)\.(?P<measure>[^.]+)$")
# Precursor quantity of a Spectronaut peptide pivot, in order of preference.
_SPECTRONAUT_PRECURSOR_QUANTITIES = ("EG.TotalQuantity (Settings)", "FG.Quantity", "FG.MS2Quantity",
                                     "PEP.Quantity", "PEP.MS2Quantity")


def is_spectronaut_long_report(columns) -> bool:
    """A Spectronaut long ("normal") report carrying protein and precursor quantities."""
    names = {str(column).strip() for column in columns}
    return ({"R.FileName", "PG.Quantity"}.issubset(names) and bool(names & {"PG.ProteinGroups", "PG.ProteinAccessions"})
            and bool(names & {"EG.ModifiedSequence", "EG.PrecursorId", "PEP.StrippedSequence"}))


def holds_protein_and_peptide_level(columns) -> bool:
    """One file the app reads as the protein table and as the peptide report."""
    return is_spectronaut_long_report(columns)


def _header_columns(path) -> list[str] | None:
    """The column names of the file at `path`, read from its header alone; None when it cannot be read."""
    try:
        path = os.fspath(path)
        extension = get_extension(path)
        if extension == "parquet":
            columns = pq.ParquetFile(path).schema_arrow.names
        elif extension == "xlsx":
            columns = pd.read_excel(path, nrows=0).columns
        else:
            columns = pd.read_csv(path, sep="," if extension == "csv" else "\t", nrows=0).columns
    except Exception:  # noqa: BLE001 - an unreadable file is reported where it is loaded
        return None
    return [str(column).strip() for column in columns]


def file_holds_protein_and_peptide_level(path) -> bool:
    """Whether the file at `path` is such a report, judged by its header alone."""
    columns = _header_columns(path)
    return columns is not None and holds_protein_and_peptide_level(columns)


def _drop_maxquant_decoys(frame: pd.DataFrame) -> pd.DataFrame:
    """MaxQuant marks reversed-sequence hits with "+" in Reverse; they are not measurements.

    Protein groups with "+" in "Only identified by site" are identified by
    modified peptides alone and are removed with them, as is usual for
    proteinGroups.txt (Tyanova et al. 2016). Contaminants ("Potential
    contaminant") stay: removing them is the user's choice in the contaminant
    filter.
    """
    def marked(column):
        if column not in frame.columns:
            return np.zeros(len(frame), dtype=bool)
        return frame[column].astype("string").str.strip().eq("+").fillna(False).to_numpy(dtype=bool)

    decoy, site_only = marked("Reverse"), marked("Only identified by site")
    site_only &= ~decoy
    if not decoy.any() and not site_only.any():
        return frame
    kept = frame.loc[~(decoy | site_only)].reset_index(drop=True)
    if decoy.any():
        kept.attrs["decoy_rows_removed"] = int(decoy.sum())
    if site_only.any():
        kept.attrs["loading_notes"] = [
            f"{int(site_only.sum()):,} protein groups identified only by a modification site "
            "(\"+\" in MaxQuant's Only identified by site column) were removed when the file was loaded."
        ]
    return kept


def _spectronaut_sequences(frame: pd.DataFrame) -> pd.DataFrame:
    """Stripped and modified sequence from EG.PrecursorId ("_PEPT[Phospho (STY)]IDE_.2") when the report lacks them."""
    source = next((column for column in ("EG.PrecursorId", "EG.ModifiedSequence", "EG.ModifiedPeptide") if column in frame.columns), None)
    if source is None:
        return frame
    text = frame[source].astype("string")
    modified = text.str.replace(r"\.\d+$", "", regex=True)
    additions = {}
    if not {"EG.ModifiedSequence", "Modified.Sequence"} & set(frame.columns):
        additions["Modified.Sequence"] = modified
    if not {"PEP.StrippedSequence", "Stripped.Sequence"} & set(frame.columns):
        additions["Stripped.Sequence"] = modified.str.replace(r"\[[^\]]*\]|_", "", regex=True)
    return frame.assign(**additions) if additions else frame


def _spectronaut_precursor_quantity(frame: pd.DataFrame) -> pd.DataFrame:
    """The precursor quantity of a long report under the name the peptide pages read."""
    if "Precursor.Quantity" in frame.columns or "FG.Quantity" in frame.columns:
        return frame
    for column in _SPECTRONAUT_PRECURSOR_QUANTITIES:
        if column in frame.columns:
            return frame.rename(columns={column: "Precursor.Quantity"})
    return frame


def _spectronaut_peptide_pivot(frame: pd.DataFrame) -> pd.DataFrame:
    """A Spectronaut peptide or precursor pivot report as a long report.

    One row per precursor and run in which it has a quantity; the run is named
    by its column stem ("[1] GFP DMSO1"). Apex and predicted retention time
    come along when the pivot carries them. A pivot without a peptide-level
    quantity is not a peptide report and is returned unchanged.
    """
    per_run: dict[str, dict[str, str]] = {}
    for column in frame.columns:
        match = _SPECTRONAUT_RUN_COLUMN.match(str(column))
        if match is not None:
            per_run.setdefault(match["run"], {})[f"{match['level']}.{match['measure']}"] = str(column)
    measures = {measure for columns in per_run.values() for measure in columns}
    quantity = next((measure for measure in _SPECTRONAUT_PRECURSOR_QUANTITIES if measure in measures), None)
    if quantity is None or any(str(column).strip() in {"R.FileName", "File.Name", "Run"} for column in frame.columns):
        return frame
    carried = {quantity: "Precursor.Quantity", "EG.ApexRT": "RT", "EG.RTPredicted": "Predicted.RT"}
    run_columns = {column for columns in per_run.values() for column in columns.values()}
    identifiers = frame[[column for column in frame.columns if str(column) not in run_columns]]
    parts = []
    for run, columns in per_run.items():
        if quantity not in columns:
            continue
        values = pd.to_numeric(frame[columns[quantity]], errors="coerce")
        measured = values.notna().to_numpy(dtype=bool)
        part = identifiers.loc[measured].copy()
        part.insert(0, "File.Name", run)
        for measure, name in carried.items():
            if measure in columns:
                part[name] = pd.to_numeric(frame.loc[measured, columns[measure]], errors="coerce").to_numpy()
        parts.append(part)
    if not parts:
        return frame
    return pd.concat(parts, ignore_index=True)


# "Control_1 Intensity", "Control_1 MaxLFQ Intensity", "Control_1 Match Type", "Control_1 Apex Retention Time"
_FRAGPIPE_RUN_COLUMN = re.compile(
    r"^(?P<run>.+?) (?P<measure>MaxLFQ Intensity|Intensity|Match Type|Apex Retention Time|Spectral Count|Localization)$")
# Mass tags of FragPipe's "Modified Sequence" ("C[57.0215]", "n[42.0106]") and
# the names the modification plot counts.
_FRAGPIPE_MODIFICATIONS = {"57.02": "Carbamidomethyl", "15.99": "Oxidation", "42.01": "Acetyl", "79.97": "Phospho"}


def is_fragpipe_peptide_table(columns) -> bool:
    """FragPipe combined_ion.tsv, combined_modified_peptide.tsv or combined_peptide.tsv (IonQuant)."""
    names = [str(column).strip() for column in columns]
    if "Peptide Sequence" not in names or {"Spectrum", "Spectrum File", "File.Name", "Run"} & set(names):
        return False
    return any(match is not None and match["measure"] == "Intensity" for match in map(_FRAGPIPE_RUN_COLUMN.match, names))


def _fragpipe_protein_names(frame: pd.DataFrame) -> pd.Series | None:
    """The entry name ("PRDX6_HUMAN") of FragPipe's "Protein" header ("sp|P30041|PRDX6_HUMAN").

    A contaminant keeps its whole header ("contam_sp|P00761|TRYP_PIG"): the
    contaminant filter reads the marker and both identifiers from it. Rows
    without an entry name (a FASTA that is not in UniProt format) keep the
    header too.
    """
    if not {"Protein", "Entry Name"}.issubset(frame.columns):
        return None
    header = frame["Protein"].astype("string").str.strip()
    entry = frame["Entry Name"].astype("string").str.strip()
    return entry.where(entry.notna() & entry.ne("") & ~header.fillna("").str.startswith("contam_"), header)


def _fragpipe_protein_table(frame: pd.DataFrame) -> pd.DataFrame:
    """FragPipe combined_protein.tsv with the protein group as ProteinNames and ProteinAccessions.

    The group is the row's protein followed by its "Indistinguishable
    Proteins", separated by semicolons like the groups of other tools, so the
    contaminant filter judges every member. Members are named like the lead
    protein (see `_fragpipe_protein_names`). Isoforms keep their accession
    ("P05408-2") and share the entry name, as in the other tools' tables.
    Nothing else changes.
    """
    names = _fragpipe_protein_names(frame)
    if names is None or "ProteinNames" in frame.columns or "ProteinAccessions" in frame.columns:
        return frame

    def members(value):
        if value is None or pd.isna(value) or not str(value).strip():
            return [], []
        entry_names, accessions = [], []
        for header in (part.strip() for part in str(value).split(",")):
            parts = header.removeprefix("contam_").split(" ", 1)[0].split("|")
            uniprot = len(parts) == 3 and parts[0] in {"sp", "tr"}
            entry_names.append(parts[2] if uniprot and not header.startswith("contam_") else header)
            accessions.append(parts[1] if uniprot else header)
        return entry_names, accessions

    others = (frame["Indistinguishable Proteins"].map(members) if "Indistinguishable Proteins" in frame.columns
              else pd.Series([([], [])] * len(frame), index=frame.index))
    accession = frame["Protein ID"].astype("string").str.strip()
    adapted = frame.copy()
    adapted.insert(0, "ProteinAccessions", [";".join([str(lead), *other[1]]) if not pd.isna(lead) else lead
                                            for lead, other in zip(accession, others)])
    # Isoforms share an entry name; a name is listed once per group, as in Spectronaut's PG.ProteinNames.
    adapted.insert(0, "ProteinNames", [";".join(dict.fromkeys([str(lead), *other[0]])) if not pd.isna(lead) else lead
                                       for lead, other in zip(names, others)])
    return adapted


def _fragpipe_peptide_table(frame: pd.DataFrame) -> pd.DataFrame:
    """A FragPipe combined ion, modified-peptide or peptide table as a long report.

    One row per table row and run with an intensity above 0 (IonQuant writes 0
    where it has no quantity); the run is the experiment name of the column
    ("Control_1"), which is also the stem of the protein table's
    "Control_1 MaxLFQ Intensity". The quantity is the "Intensity" column, in
    all three tables, including values transferred by match between runs; the
    match type comes along as Match.Type. The apex retention time (seconds in
    the ion table) is given in minutes. FragPipe has no predicted retention
    time, so the retention-time plot says that the report lacks it. Besides
    "MS/MS" and "MBR" a newer IonQuant writes "MBC": an ion quantified
    without a spectrum of its own in a run where the same modified peptide
    was identified at another charge (all 6,235 such values of PXD082879).
    It is a quantity like the others. The protein is the one FragPipe
    assigned the peptide to, named by its entry name.
    """
    names = {column: str(column).strip() for column in frame.columns}
    per_run: dict[str, dict[str, object]] = {}
    for column, name in names.items():
        match = _FRAGPIPE_RUN_COLUMN.match(name)
        if match is not None:
            per_run.setdefault(match["run"], {})[match["measure"]] = column
    per_run = {run: columns for run, columns in per_run.items() if "Intensity" in columns}
    run_columns = {column for columns in per_run.values() for column in columns.values()}
    identifiers = frame[[column for column in frame.columns if column not in run_columns]].rename(columns=names)
    identifiers = identifiers.rename(columns={"Peptide Sequence": "Stripped.Sequence", "Charge": "Precursor.Charge"}).copy()
    protein_names = _fragpipe_protein_names(identifiers)
    if protein_names is not None:
        identifiers.insert(0, "ProteinNames", protein_names)
    if "Modified Sequence" in identifiers.columns:
        def named(tag):
            return "[" + _FRAGPIPE_MODIFICATIONS.get(f"{float(tag[1]):.2f}", tag[1]) + "]"
        identifiers["Modified.Sequence"] = identifiers["Modified Sequence"].astype("string").str.replace(
            r"\[(\d+(?:\.\d+)?)\]", named, regex=True)
        if "Precursor.Charge" in identifiers.columns:
            charge = pd.to_numeric(identifiers["Precursor.Charge"], errors="coerce").astype("Int64").astype("string").fillna("")
            identifiers["Precursor.Id"] = identifiers["Modified.Sequence"].astype("string") + charge
    parts = []
    for run, columns in per_run.items():
        values = pd.to_numeric(frame[columns["Intensity"]], errors="coerce")
        measured = (values.notna() & values.gt(0)).to_numpy(dtype=bool)
        part = identifiers.loc[measured].copy()
        part.insert(0, "File.Name", run)
        part["Precursor.Quantity"] = values.to_numpy()[measured]
        if "Apex Retention Time" in columns:
            part["RT"] = pd.to_numeric(frame.loc[measured, columns["Apex Retention Time"]], errors="coerce").to_numpy() / 60.0
        if "Match Type" in columns:
            part["Match.Type"] = frame.loc[measured, columns["Match Type"]].to_numpy()
        parts.append(part)
    return pd.concat(parts, ignore_index=True)


def _maxquant_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    """MaxQuant evidence.txt under the peptide report's column names.

    The run is the Experiment when the file has one (that is the name
    proteinGroups.txt uses in "LFQ intensity <experiment>"), otherwise the raw
    file. MaxQuant has no predicted retention time, so the retention-time plot
    stays unavailable.
    """
    renamed = frame.rename(columns={
        "Sequence": "Stripped.Sequence", "Modified sequence": "Modified.Sequence", "Intensity": "Precursor.Quantity",
        "Retention time": "RT", "Charge": "Precursor.Charge", "Leading razor protein": "ProteinNames"})
    run = "Experiment" if "Experiment" in renamed.columns else "Raw file"
    renamed.insert(0, "File.Name", renamed[run].astype(str))
    if "Precursor.Charge" in renamed.columns:
        charge = pd.to_numeric(renamed["Precursor.Charge"], errors="coerce").astype("Int64").astype("string").fillna("")
        renamed["Precursor.Id"] = renamed["Modified.Sequence"].astype(str) + charge
    return renamed


def _spectronaut_site_report(frame: pd.DataFrame) -> pd.DataFrame:
    """Spectronaut's own PTM site report under the app's phosphosite columns.

    Spectronaut collapses precursors to sites itself (PTM.CollapseKey,
    PTM.Quantity). The app's phosphosite pages read the columns of its own
    collapse: PTM_Collapse_key ("<protein>~<gene>_<residue><position>_M<n>"),
    Protein_group, Gene_group, UPD_seq and PTM_localization. Only phospho rows
    are kept. UPD_seq is the flanking window with the site in lower case and
    marked "*"; the gene is the first PG.Genes seen for the protein, or the
    accession when the report has no gene column. PTM_localization is the
    highest PTM.SiteProbability of the site, and missing when the report was
    exported without it. A long report (one row per site and run) is pivoted
    to one column per run.
    """
    data = frame
    if "PTM.ModificationTitle" in data.columns:
        phospho = data["PTM.ModificationTitle"].astype("string").str.startswith("Phospho").fillna(False).to_numpy(dtype=bool)
        if not phospho.any():
            raise ValueError("The Spectronaut PTM report contains no phosphorylation sites. Export a report containing Phospho modifications.")
        data = data.loc[phospho]
    per_run = {str(column): match for column in data.columns
               if (match := _SPECTRONAUT_RUN_COLUMN.match(str(column))) is not None and match["level"] == "PTM"}
    probability = None
    if "R.FileName" in data.columns and "PTM.Quantity" in data.columns:
        runs = list(dict.fromkeys(data["R.FileName"].astype(str)))
        values = (data.assign(**{"PTM.Quantity": pd.to_numeric(data["PTM.Quantity"], errors="coerce")})
                  .pivot_table(index="PTM.CollapseKey", columns="R.FileName", values="PTM.Quantity", aggfunc="first",
                               sort=False, dropna=False).reindex(columns=runs))
        values.columns = [f"{run}.PTM.Quantity" for run in values.columns]
        if "PTM.SiteProbability" in data.columns:
            probability = pd.to_numeric(data["PTM.SiteProbability"], errors="coerce").groupby(
                data["PTM.CollapseKey"].to_numpy(), sort=False).max()
        described = [column for column in data.columns if str(column).startswith(("PTM.", "PG.")) and str(column) not in
                     {"PTM.Quantity", "PTM.SiteProbability", "PG.Quantity", "PG.Qvalue", "PG.QValue (Run-Wise)"}]
        data = data[described].drop_duplicates("PTM.CollapseKey").set_index("PTM.CollapseKey", drop=False)
        data = data.join(values).reset_index(drop=True)
        if probability is not None:
            probability = data["PTM.CollapseKey"].map(probability)
    else:
        chances = [column for column, match in per_run.items() if match["measure"] == "SiteProbability"]
        if chances:
            probability = data[chances].apply(pd.to_numeric, errors="coerce").max(axis=1)
        data = data.drop(columns=[column for column, match in per_run.items() if match["measure"] != "Quantity"])
    data = data.reset_index(drop=True)
    protein = data["PTM.ProteinId"].astype(str).str.strip()
    if "PG.Genes" in data.columns:
        genes = data["PG.Genes"].astype("string").str.split(";").str[0].str.strip()
        gene = genes.groupby(protein.to_numpy(), sort=False).transform("first").fillna(protein).where(lambda value: value.ne(""), protein)
    else:
        gene = protein
    gene = gene.astype(str).str.replace("_", "#", regex=False)
    residue = data["PTM.SiteAA"].astype(str).str.strip()
    position = pd.to_numeric(data["PTM.SiteLocation"], errors="coerce").astype("Int64").astype("string")
    multiplicity = pd.to_numeric(data["PTM.Multiplicity"], errors="coerce").astype("Int64").astype("string") \
        if "PTM.Multiplicity" in data.columns else pd.Series("1", index=data.index, dtype="string")
    window = data["PTM.FlankingRegion"].astype(str) if "PTM.FlankingRegion" in data.columns else pd.Series("", index=data.index)
    centre = window.str.len() // 2
    site = [text[:middle] + text[middle].lower() + "*" + text[middle + 1:] if text and text != "nan" else ""
            for text, middle in zip(window.tolist(), centre.tolist())]
    data.insert(0, "PTM_Collapse_key", protein + "~" + gene + "_" + residue + position.fillna("") + "_M" + multiplicity.fillna("1"))
    data.insert(1, "Protein_group", protein)
    data.insert(2, "Gene_group", gene)
    data.insert(3, "UPD_seq", site)
    data.insert(4, "PTM_localization", np.nan if probability is None else pd.to_numeric(probability, errors="coerce").to_numpy())
    return data[position.notna().to_numpy(dtype=bool)].reset_index(drop=True)


def adapt_report(frame: pd.DataFrame, peptide: bool = False) -> pd.DataFrame:
    """Reshape the reports of other tools the app has an adapter for (see above)."""
    columns = [str(column).strip() for column in frame.columns]
    frame = _drop_maxquant_decoys(frame)
    if peptide:
        if {"Sequence", "Modified sequence", "Raw file"}.issubset(columns):
            return _maxquant_evidence(frame)
        if is_fragpipe_peptide_table(columns):
            return _fragpipe_peptide_table(frame)
        frame = _spectronaut_peptide_pivot(frame)
        if any(column.startswith(("EG.", "PEP.")) for column in columns):
            frame = _spectronaut_precursor_quantity(_spectronaut_sequences(frame))
        return frame
    if {"PTM.CollapseKey", "PTM.ProteinId", "PTM.SiteAA", "PTM.SiteLocation"}.issubset(columns):
        return _spectronaut_site_report(frame)
    if is_spectronaut_long_report(columns):
        return _spectronaut_protein_matrix(frame)
    if {"Protein", "Protein ID", "Entry Name"}.issubset(columns) and any(column.endswith(" Intensity") for column in columns):
        return _fragpipe_protein_table(frame.rename(columns=lambda column: str(column).strip()))
    return frame


def _spectronaut_protein_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """The protein groups of a Spectronaut long report as a matrix: PG.Quantity per run.

    PG.Quantity is repeated on every precursor row of a protein group and run;
    the first value is taken. Columns are named "<run>.PG.Quantity" like the
    pivot report's.
    """
    key = "PG.ProteinGroups" if "PG.ProteinGroups" in frame.columns else "PG.ProteinAccessions"
    runs = list(dict.fromkeys(frame["R.FileName"].astype(str)))
    values = (frame.assign(**{"PG.Quantity": pd.to_numeric(frame["PG.Quantity"], errors="coerce")})
              .pivot_table(index=key, columns="R.FileName", values="PG.Quantity", aggfunc="first", sort=False, dropna=False)
              .reindex(columns=runs))
    values.columns = [f"{run}.PG.Quantity" for run in values.columns]
    described = [column for column in frame.columns if str(column).startswith("PG.") and str(column) not in
                 {"PG.Quantity", "PG.Qvalue", "PG.QValue (Run-Wise)", "PG.Cscore", "PG.Cscore (Run-Wise)", "PG.PEP",
                  "PG.PEP (Run-Wise)", "PG.Pvalue", "PG.PValue (Run-Wise)", "PG.RunEvidenceCount", "PG.IsSingleHit",
                  "PG.Coverage", "PG.MS1Quantity", "PG.MS2Quantity"}]
    groups = frame[described].drop_duplicates(key).set_index(key, drop=False)
    return groups.join(values).reset_index(drop=True)


def sanitize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = make_columns_unique([str(c) for c in df.columns])

    for col in df.columns:
        if not pd.api.types.is_numeric_dtype(df[col]):
            # Keep missing cells missing. On pandas < 3, `astype(str)` writes
            # the words "nan"/"None" into them, and those then pass as gene or
            # protein names -- so the same file carried a gene called "nan" on
            # one pandas and no gene on another.
            column = df[col]
            df[col] = column.astype(str).where(column.notna(), None)

    return df


_CHUNK_ROWS = 50_000


def _rewind(fileobj: BinaryIO) -> None:
    try:
        fileobj.seek(0)
    except Exception:  # noqa: BLE001 - an unseekable stream is simply read once
        pass


def _stream_length(fileobj: BinaryIO) -> int | None:
    """Byte length of a seekable stream, or None when it cannot be measured."""
    try:
        current = fileobj.tell()
        fileobj.seek(0, 2)
        length = fileobj.tell()
        fileobj.seek(current)
        return int(length)
    except Exception:  # noqa: BLE001 - an unseekable stream is read in chunks instead
        return None


def _frame_bytes(frame: pd.DataFrame) -> int:
    """Materialised size of a frame, following object dtypes to the real strings."""
    try:
        return int(frame.memory_usage(deep=True).sum())
    except Exception:  # noqa: BLE001 - an exotic dtype should not defeat the budget
        return int(frame.memory_usage(deep=False).sum())


_NUL_SCAN_BYTES = 64 * 1024


def _reject_embedded_nul(fileobj: BinaryIO, filename: str) -> None:
    """Refuse delimited text carrying a NUL byte.

    pandas' C parser treats NUL as a terminator, so a header written as
    "Protein\x00Names" arrives as "Protein" -- the column is silently renamed and
    every later lookup for the real name misses. A NUL in what is meant to be
    UTF-8 text almost always means the file is UTF-16 or was truncated mid-write,
    so it is better refused than quietly repaired.

    Only the first chunk is scanned: that covers the header, where this actually
    bites, without reading a file the streaming budget exists to avoid.
    """
    try:
        fileobj.seek(0)
        head = fileobj.read(_NUL_SCAN_BYTES)
    except Exception:  # noqa: BLE001 - an unseekable stream is parsed as it comes
        return
    finally:
        _rewind(fileobj)
    if bytes([0]) in head:
        raise ValueError(
            f"Could not read {filename}: the file contains a NUL byte where text was "
            "expected, which usually means it is UTF-16 or was truncated. Re-save it "
            "as UTF-8 text."
        )


def _read_delimited_within_budget(fileobj: BinaryIO, *, sep: str, filename: str, select_column=None) -> pd.DataFrame:
    """Read a delimited table, stopping as soon as it exceeds a budget.

    Three things are bounded, because a file can grow in three directions and a
    check on one of them says nothing about the others: the column count, the
    cell count, and the bytes the frame actually occupies. Reading in chunks caps
    the cost of an oversized file at one chunk beyond the limit rather than the
    whole thing.
    """
    # Delimited text is not compressed, so its length on the wire is an honest
    # floor for what the frame will occupy -- a value costs at least its
    # characters, and usually more once it becomes a float or a Python string.
    # Checking it first stops a file of few but enormous rows, which row-based
    # chunking cannot: 400 rows is one chunk, so the chunk loop would read the
    # whole thing before it got to weigh any of it.
    raw_length = _stream_length(fileobj)
    if raw_length is not None:
        ensure_frame_bytes_within_budget(raw_length, label="Upload")
    _reject_embedded_nul(fileobj, filename)

    source_width = None
    if select_column is not None:
        headers = pd.read_csv(fileobj, sep=sep, nrows=0).columns
        source_width = len(headers)
        ensure_columns_within_budget(source_width, label="Upload")
        _rewind(fileobj)
    chunks: list[pd.DataFrame] = []
    rows = 0
    columns = 0
    nbytes = 0
    for chunk in pd.read_csv(fileobj, sep=sep, chunksize=_CHUNK_ROWS, usecols=select_column):
        columns = source_width or len(chunk.columns)
        # Checked on every chunk, not just the first: the width is known from the
        # header, and a wide file puts its whole width in that first chunk.
        ensure_columns_within_budget(columns, label="Upload")
        rows += len(chunk)
        ensure_table_within_cell_budget(rows, columns, label="Upload")
        nbytes += _frame_bytes(chunk)
        ensure_frame_bytes_within_budget(nbytes, label="Upload")
        chunks.append(chunk)
    if not chunks:
        _rewind(fileobj)
        return pd.read_csv(fileobj, sep=sep, usecols=select_column)
    return pd.concat(chunks, ignore_index=True) if len(chunks) > 1 else chunks[0]


def _ensure_workbook_within_budget(fileobj: BinaryIO) -> None:
    """Check an xlsx against the budgets before openpyxl materialises anything.

    An xlsx is a zip, and a zip's central directory carries each member's
    uncompressed size, so the expanded cost is knowable without decompressing.
    That is the measure the shape-only check missed: a workbook of a few cells
    holding megabytes of text each is small by row and column count.
    """
    try:
        import zipfile

        fileobj.seek(0)
        with zipfile.ZipFile(fileobj) as archive:
            expanded = sum(int(entry.file_size) for entry in archive.infolist())
        ensure_frame_bytes_within_budget(expanded, label="Upload")
    except ValueError:
        raise
    except Exception:  # noqa: BLE001 - a malformed zip is pandas' to report
        pass
    finally:
        _rewind(fileobj)

    try:
        import openpyxl

        workbook = openpyxl.load_workbook(fileobj, read_only=True)
    except ValueError:
        raise
    except Exception:  # noqa: BLE001 - let pandas report a malformed workbook
        _rewind(fileobj)
        return
    try:
        for sheet in workbook.worksheets:
            ensure_columns_within_budget(sheet.max_column or 0, label="Upload")
            ensure_table_within_cell_budget(
                sheet.max_row or 0, sheet.max_column or 0, label="Upload"
            )
    finally:
        workbook.close()
        _rewind(fileobj)


def _read_parquet_within_budget(parquet_file, columns=None) -> pd.DataFrame:
    """Materialise a parquet file in batches, stopping if it exceeds the budget."""
    frames: list[pd.DataFrame] = []
    nbytes = 0
    for batch in parquet_file.iter_batches(batch_size=_CHUNK_ROWS, columns=columns):
        frame = batch.to_pandas()
        nbytes += _frame_bytes(frame)
        ensure_frame_bytes_within_budget(nbytes, label="Upload")
        frames.append(frame)
    if not frames:
        return parquet_file.read(columns=columns).to_pandas()
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]


def _parse_dataframe(filename: str, fileobj: BinaryIO, select_column=None, peptide: bool = False) -> pd.DataFrame:
    ext = get_extension(filename)

    if ext in {"csv", "tsv", "txt"}:
        df = _read_delimited_within_budget(fileobj, sep="," if ext == "csv" else "\t", filename=filename, select_column=select_column)
    elif ext == "xlsx":
        # xlsx is a compressed container too, so the shape is checked from the
        # worksheet dimensions before openpyxl materialises the cells.
        _ensure_workbook_within_budget(fileobj)
        df = pd.read_excel(fileobj, usecols=select_column)
    elif ext == "parquet":
        # Parquet carries its own shape, so the budget is checked before the frame
        # is materialised. Doing it afterwards would mean allocating the very
        # thing the limit exists to prevent.
        parquet_file = pq.ParquetFile(fileobj)
        metadata = parquet_file.metadata
        ensure_columns_within_budget(metadata.num_columns, label="Upload")
        ensure_table_within_cell_budget(
            metadata.num_rows, metadata.num_columns, label="Upload"
        )
        # The footer records each row group's uncompressed size. That is an
        # honest measure when values are distinct, and an undercount when they
        # repeat, because parquet dictionary encodes them -- a file reporting
        # 20 MB can still build 600 MB of separate Python strings. So it is a
        # cheap first refusal, not the only one.
        expanded = sum(
            metadata.row_group(index).total_byte_size
            for index in range(metadata.num_row_groups)
        )
        ensure_frame_bytes_within_budget(expanded, label="Upload")
        # The footer's uncompressed size is a floor, not the cost: identical
        # values are dictionary encoded, so a file that reports 20 MB can expand
        # past 600 MB once every row gets its own Python string. Read in batches
        # and measure what is actually built, stopping one batch past the limit.
        columns = [c for c in parquet_file.schema_arrow.names if select_column(c)] if select_column else None
        df = _read_parquet_within_budget(parquet_file, columns)
    else:
        raise ValueError(f"Unsupported file type: {ext}")

    ensure_columns_within_budget(len(df.columns), label="Upload")
    ensure_table_within_cell_budget(len(df), len(df.columns), label="Upload")
    ensure_frame_bytes_within_budget(_frame_bytes(df), label="Upload")

    df = adapt_report(df, peptide)
    # A pivot turned into a long report grows by its number of runs.
    ensure_frame_bytes_within_budget(_frame_bytes(df), label="Upload")
    removed, notes = df.attrs.get("decoy_rows_removed"), df.attrs.get("loading_notes")
    df = rename_cols(df)
    df = sanitize_dataframe(df)
    if removed:
        df.attrs["decoy_rows_removed"] = removed
    if notes:
        df.attrs["loading_notes"] = list(notes)
    return df


_PARSE_HINTS = (
    ("no columns to parse", "the file is empty or has no header row"),
    ("codec can't decode", "the file is not valid UTF-8 text -- if it was exported from "
                           "Excel or instrument software, re-save it as UTF-8"),
    ("error tokenizing", "its rows do not all have the same number of columns"),
    ("excel file format cannot be determined", "it is named .xlsx but is not a valid workbook"),
)


def _readable_parse_error(filename: str, error: Exception) -> str:
    """Turn a parser's own words into something a biologist can act on."""
    text = str(error)
    lowered = text.lower()
    for needle, explanation in _PARSE_HINTS:
        if needle in lowered:
            return f"Could not read {filename}: {explanation}."
    return f"Could not read {filename}: {text}"


def read_dataframe(filename: str, fileobj: BinaryIO, *, select_column: Callable[[str], bool] | None = None,
                   peptide: bool = False) -> pd.DataFrame:
    """Read an uploaded table, refusing one that is not a table at all.

    Two things go wrong here that used to be reported badly or not at all. A file
    that parses but holds no measurements -- a saved error page, a truncated
    export, a header with no rows -- was accepted with a success preview, and the
    user found out several screens later. And a parser's own exception text was
    passed straight through, so "Error tokenizing data. C error: Expected 3
    fields in line 3, saw 5" was the whole explanation offered.
    """
    try:
        frame = _parse_dataframe(filename, fileobj, select_column, peptide)
    except ValueError as error:
        # The upload budgets raise ValueError too, and those messages are already
        # written for the reader, so they must pass through untouched.
        if str(error).startswith(("Upload ", "Unsupported file type", "Could not read")):
            raise
        raise ValueError(_readable_parse_error(filename, error)) from error
    except Exception as error:  # noqa: BLE001 - a parser failure is the user's problem to fix
        raise ValueError(_readable_parse_error(filename, error)) from error

    if len(frame) == 0:
        raise ValueError(
            f"{filename} has no data rows. Check the export finished and that the "
            "header is the first line of the file."
        )
    if len(frame.columns) < 2:
        raise ValueError(
            f"{filename} has only {len(frame.columns)} column. A report table needs "
            "an identifier column and at least one measurement column -- check the "
            "file is not a saved web page and that the separator is correct."
        )
    return frame
