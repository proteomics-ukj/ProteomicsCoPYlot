# Changelog

## 1.0.0-beta.2 — unreleased

- Export and restore all original metadata profiles in one Excel workbook,
  preserving profile names, sample grouping, and pairing.
- Improve report editing, add custom sample labels, and align phosphosite
  coverage sections with their plot settings.
- Offer report notebooks only when every displayed section was captured;
  allow large exports beyond the soft notebook cache targets.
- Reuse peptide preparation and ANOVA calculations when restyling figures,
  and space dense heatmap row labels without removing data rows.
- Use static scatter plots by default for peptide retention-time comparisons
  and remove the interactive retention-time option from the graphical interface.
- Simplify plot and filtering controls and standardize default random seeds
  to 1 across the application and notebook helpers. Analyses using defaults
  can therefore differ from Beta 1 when they involve random sampling or imputation.
- Update Python, backend, frontend, citation, and API versions to Beta 2.

## 1.0.0-beta.1 — unreleased

- Mark the first public release as Beta 1 across Python, backend, frontend,
  citation, and API metadata.
- Use `1.0.0b1` for Python packages and `1.0.0-beta.1` for frontend and release
  directory names; the GitHub tag is `v1.0.0-beta.1`.
- Build distribution files from the merged `main` commit before publishing.
