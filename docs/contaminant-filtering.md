# Contaminant filtering

In **Filter Settings**, enable **Filter contaminants**, select categories, then press
the existing **Apply filter** button. It is off by default. Enabling it initially
selects all five categories from the bundled Cambridge `crap_ccp_metadata.csv` plus
**Upstream-flagged contaminants**; any nonempty combination can be selected. This sixth category
removes rows marked `Potential contaminant = +` by MaxQuant, independently of cRAP membership
and accession matching. It also removes groups whose supplied protein-name/accession members
are all prefixed `CON__` (MaxQuant), `contam_` (FragPipe/Philosopher) or `cRAP-` (a cRAP FASTA
searched alongside the proteome, as in DIA-NN `cRAP-TRYP_PIG`). Within one identifier column
every member must carry the prefix: a group that pairs a prefixed entry with its unprefixed twin
(`K2C1_HUMAN;cRAP-K2C1_HUMAN`) is left to the database categories. A column without any prefix,
such as a bare accession column beside prefixed names, may name the same proteins. These prefixes are also
recognized when retained in DIA-NN or Spectronaut exports. No contamination is inferred
from a gene name or an arbitrary substring. Mixed or incomplete prefixed groups stay unless the
explicit row flag is set. `Reverse` and `Only identified by site` are not contaminant flags;
rows carrying them are removed when a MaxQuant table is loaded.
Category 6 can be selected alone or together with database categories; existing selections
are preserved; new default configurations select all six categories. The 14 records classified by scientific
decision participate in their assigned categories just like the 111 GPM matches.
If a metadata record has several category IDs separated by semicolons, selecting
any one of those categories includes that protein.

The valid-value filter runs first, followed by contaminant removal. Normalization
is applied separately and excludes contaminant rows from its fit, including when
fitting on complete features. Original uploads and unfiltered matrices remain
available; use the **filtered** source for analyses that should exclude contaminants.
If every row is excluded, the filtered table is empty and normalization is refused.

## Matching

- Exact, case-sensitive protein entry names (such as `ALBU_HUMAN`) and complete
  protein names are matched. Gene symbols, substrings and sequence similarity are
  not used. Whitespace around identifiers is ignored.
- **Advanced matching → Also match exact protein accessions when available** is
  enabled by default. It supports accession columns and accession identifiers in
  the protein-name column, including UniProt `sp|accession|entry` / `tr|...` headers.
  The original Cambridge `sp|cRAPnnn|accession|entry` headers are also recognised.
  Disable it to use names only. Isoform suffixes are not stripped. The placeholder
  accession `000000` never matches; those two records require their entry names.
- MaxQuant `CON__`, FragPipe `contam_` and cRAP-FASTA `cRAP-` prefixes are stripped from each member
  only for matching. For example, `CON__P02769` and
  `contam_sp|P02769|ALBU_BOVIN` can match the selected reference entries, as can
  `cRAP-ALBU_BOVIN` / `cRAP-P02769`; a pair such as `K2C1_HUMAN;cRAP-K2C1_HUMAN` is one
  fully matched group.
  Decoy prefixes are not stripped. Prefixes and `Potential contaminant = +`
  alone never remove a row through categories 1–5; a selected-category database match is still required.
  Category 6 uses the independent upstream-flag rules described above.
- Protein groups use semicolon-separated members. Every member must match the
  union of selected categories; a contaminant plus an unmatched protein stays.
  Missing or empty group members do not match.
- When both names and accessions are supplied, the complete accession group must
  match. This can resolve unfamiliar display names, but a shorter accession list
  cannot remove a larger, partly unmatched name group. Conflicting or incomplete
  accession groups therefore retain rows. Without accessions, names are used.
- Canonical name columns are `ProteinNames` or `Phosphoprotein`; accession columns
  are `ProteinAccessions` or `Protein_group`. Upload aliases are mapped as usual.
  Selecting database categories without a usable identifier column produces an error.
  MaxQuant contaminant flags alone do not trigger this database filter.
- Decoy rows of a MaxQuant table (`Reverse = +`) and protein groups identified only by a
  modification site (`Only identified by site = +`) are not contaminants and are not
  this filter's business: they are removed when the file is loaded, and the
  annotation states how many.

Peptide plots and their notebook exports use the unfiltered peptide report.
Protein and phosphosite filter settings are independent. Phosphoprotein aggregation inherits the upstream phosphosite
filter, as it does the existing valid-value filter.

The filter records its post-valid-value removal count, selected categories,
matching mode, database version and metadata SHA-256 in the annotation messages.
Settings are saved with the annotation/session. Notebook exports embed the exact
reference lookup, attribution/licence notices and matching code so the original upload can be replayed without
downloading a database.

Shared viewer configurations accept `contaminants_protein` and
`contaminants_phospho` objects, for example
`{"enabled": true, "categories": [1, 2], "matchAccessions": true}`. Omitted settings
default to off, including in older saved sessions.

## Standalone Python and R examples

[`examples/example.ipynb`](../examples/example.ipynb) and
[`examples/example.R`](../examples/example.R) show optional filtering between
`filter_data` and `normalize_data`. Both examples enable `filter_contaminants`
with all six categories selected. The helper itself defaults to off, with selectable
category IDs and optional accession matching. The helper
preserves Python indices / R row names, so normalization can fit exactly the
surviving rows. The examples stop if no rows remain for fitting. They also include
an enabled mixed-group demonstration and independent phosphosite settings.

A listed protein can be a genuine component of a sample. Select categories that
fit the experiment; a database match alone is not proof of contamination.
For sources, attribution, category provenance and licence notices, see the
[database documentation](../data/db/README.md) and
[Cambridge notice](../data/db/CAMBRIDGE_NOTICE.md).
