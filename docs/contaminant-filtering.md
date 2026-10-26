# Contaminant filtering

In **Filter Settings**, enable **Filter contaminants**, select categories, then press
the existing **Apply filter** button. It is off by default. Enabling it initially
selects all five categories from the bundled Cambridge `crap_ccp_metadata.csv`;
any nonempty combination can be selected. The 14 records classified by scientific
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
- MaxQuant `CON__` and FragPipe `contam_` prefixes are stripped from each member
  only for matching. For example, `CON__P02769` and
  `contam_sp|P02769|ALBU_BOVIN` can match the selected reference entries.
  Decoy prefixes are not stripped. Prefixes and `Potential contaminant = +`
  alone never remove a row; a selected-category database match is still required.
- Protein groups use semicolon-separated members. Every member must match the
  union of selected categories; a contaminant plus an unmatched protein stays.
  Missing or empty group members do not match.
- When both names and accessions are supplied, the complete accession group must
  match. This can resolve unfamiliar display names, but a shorter accession list
  cannot remove a larger, partly unmatched name group. Conflicting or incomplete
  accession groups therefore retain rows. Without accessions, names are used.
- Canonical name columns are `ProteinNames` or `Phosphoprotein`; accession columns
  are `ProteinAccessions` or `Protein_group`. Upload aliases are mapped as usual.
  Enabling the filter without a usable identifier column produces an error.
  MaxQuant contaminant/decoy flags alone do not trigger this database filter.

Protein settings also apply to peptide rows through their parent-protein
annotations, using the same category and matching rules. They do not inherit the
protein valid-value filter. Settings follow the active metadata profile; without
a protein annotation, peptide contaminant filtering is off. Phosphosite settings
are independent. Phosphoprotein aggregation inherits the upstream phosphosite
filter, as it does the existing valid-value filter.

The filter records its post-valid-value removal count, selected categories,
matching mode, database version and metadata SHA-256 in the annotation messages.
Settings are saved with the annotation/session. Notebook exports embed the exact
reference lookup, attribution/licence notices and matching code so the original upload can be replayed without
downloading a database. Peptide exports capture the inherited protein settings.

Shared viewer configurations accept `contaminants_protein` and
`contaminants_phospho` objects, for example
`{"enabled": true, "categories": [1, 2], "matchAccessions": true}`. Omitted settings
default to off, including in older saved sessions.

## Standalone Python and R examples

[`examples/example.ipynb`](../examples/example.ipynb) and
[`examples/example.R`](../examples/example.R) show optional filtering between
`filter_data` and `normalize_data`. Both expose `filter_contaminants`, defaulting
to off, with selectable category IDs and optional accession matching. The helper
preserves Python indices / R row names, so normalization can fit exactly the
surviving rows. The examples stop if no rows remain for fitting. They also include
an enabled mixed-group demonstration and independent phosphosite settings.

A listed protein can be a genuine component of a sample. Select categories that
fit the experiment; a database match alone is not proof of contamination.
For sources, attribution, category provenance and licence notices, see the
[database documentation](../data/db/README.md) and
[Cambridge notice](../data/db/CAMBRIDGE_NOTICE.md).
