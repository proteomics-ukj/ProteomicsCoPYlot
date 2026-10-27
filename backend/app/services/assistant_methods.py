"""Read-only method descriptions; never inspect or mutate an uploaded dataset."""
from app.schemas.assistant import Evidence

# "What is the best order to go through the analysis?" (2026-10-09) had only page-level
# method texts to draw on and was answered with volcano settings alone.
WORKFLOW_GUIDE = "A typical path through Proteomics CoPYlot, as a suggestion; every setting stays the user's choice. 1. Data / Upload and Data / Annotation: load the matrix, assign samples to conditions, declare whether values are log2, and apply the valid-value filter followed by optional Cambridge cRAP contaminant filtering (off by default; selected categories, exact protein names or accessions, every group member must match; no gene matching and no removal based on search-engine flags alone). Recognised CON__ and contam_ prefixes are stripped only for database matching. Protein settings also govern peptides; phosphosite settings are independent. Then (as a second step) a normalization fitted on the rows the filter kept; none are imposed. 2. QC Pipeline: Coverage Plot (identifications per sample), Histogram Intensity and Boxplot Intensity (are samples comparable in level and spread), CV Plot (replicate precision), Principal Component Analysis and Correlation Plot (do replicates group by condition; does a sample stand apart). 3. Completeness: Missing Value Plot and Heatmap, to see whether missingness differs between conditions before deciding whether to impute (Data / Imputation is optional). 4. Statistical Analysis: Volcano Plot for two conditions (choose the contrast, test and thresholds), ANOVA for more than two conditions (these are the tool's designs: two-group tests, control-corrected contrasts and a one-way ANOVA, with pairing as the only blocking; batches, covariates and repeated measures are not modelled), GO Enrichment (ORA) and Pathway Heatmap to interpret hit lists, Replicate & Variance Projection as a what-if view of the same comparison at other replicate counts (not a power calculation). 5. Single Protein for individual proteins or sites; Phospho-specific pages for phosphosite data (localization, STY, Kinase Library inputs); Comparison for selected samples or conditions within the chosen dataset level; Summary for tables, text and the report. Downloads are the user's own (the assistant never presses them): every plot has a format choice and Download Plot; result and summary tables (including Comparison Venn and Pearson, phosphosite summaries, Volcano, ANOVA and GO Enrichment) have Download CSV above the table, which Excel opens. For a request for both a figure and its numbers, explain the plot download and the table CSV separately; notebook export is the reproducibility option, not a substitute for explaining the CSV; Export notebook appears beside the download controls of supported plot or table outputs in analysis mode; it is disabled until an exportable output is rendered, and is not available in viewer mode. It downloads a ZIP containing a Python notebook, requirements.txt, a manifest, the captured input tables and resources needed to reproduce that output, settings and the original output for comparison. Extract and share the whole ZIP, not the notebook alone. The recipient also needs the original uploaded data file(s): set data_path as instructed in the notebook. Captured input matrices are verification references, not substitutes for the original upload. Use the Python version required by the exported notebook, install requirements.txt and run the cells in order. The archive contains analysis data, so inspect it before sharing. Summary / Report builds the report."

METHODS = {
    "data.translator": "ID Translator maps identifiers in a loaded protein or phospho table. Translate and Preview appends mapped values as new columns and writes the result back to the active dataset for downstream annotation; it is not a detached preview. Choose the input column and source/target identifier databases manually. The assistant can open and explain this page but cannot run translation or overwrite data. Download Translated Table saves the result manually.",
    "data.conditionPalette": "Colors defines the saved condition palette for the selected dataset level; plots use that palette. Changing a colour changes appearance, not measurements. On this page the colour fields, reset and copy actions are manual-only: choose the dataset level, edit a condition swatch or its hex field, and complete a valid colour. The assistant can open and explain this page but cannot apply its palette changes. A request to preserve measurements does not prevent styling elsewhere.",
    "external.gonetView": "GONetView visualizes Gene Ontology term networks and gene connections to biological processes, molecular functions and cellular components. It is not a STRING protein-interaction view. The CoPYlot page contains an overview and links; opening that page does not upload the active protein table. Its overview image is fetched from GitHub, so opening the page can make an external network request. Open GONetView launches the separate website; inspect that tool before choosing what to share. The assistant only navigates to this information page; it does not upload data or operate the external website.",
    "external.peptideCollapse": "Peptide Collapse accepts a manually chosen Spectronaut TSV and localization cutoff, collapses phosphosite profiles and downloads the result. These are processing choices: the assistant can navigate here and explain them but cannot choose a file, change the cutoff, run collapse or download on your behalf. Select the file and settings yourself, then use Run Peptide Collapse.",
    "data.distribution": "Pools sample values from the distribution service's selected source and reports descriptive statistics, with population standard deviation. A QQ plot compares a distribution with a normal reference; numerical summaries alone do not establish normality, batch effects or an imputation mechanism. A non-normal distribution or QQ departure alone does not show that normalization is needed: normalization concerns comparability across samples, not forcing a pooled distribution to be normal. Review sample distributions and the experimental design; processing remains your choice.",
    "stats.simulation": "A what-if projection, not a power calculation: the observed mean fold changes are held fixed and the independent-sample t-test is recomputed for an assumed variance factor and number of replicates per group. Nothing is drawn at random. The observed fold changes include noise, so the projection is optimistic for effects that look large, and it does not estimate the probability of detecting a true effect. Interpret the selected pooled-variance/Welch test and raw/BH-adjusted p-values explicitly.",
    "phospho.ksea": "Prepares testable phosphosite fold changes and p-values for an external Kinase Library tool. These are inputs, not calculated kinase activity. Uploaded external kinase plots have separate provenance and are not evidence of activity from the export table; substrate counts and kinase-level scores require the external result data.",
    "analysis.canvas": "The canvas composes multiple analysis tiles and reusable entity selections. It is not a single statistical test. Interpret each tile using its own source, contrast and method; the composition itself does not supply a combined scientific conclusion.",
    "stats.anova": "One-way independent-groups ANOVA tests whether all condition means are equal, assuming equal population variances. Rows need sufficient replicates in every condition; undefined tests are excluded. The selected metric is raw p or BH-adjusted p across finite tests. Omnibus significance does not identify a particular condition pair or direction. Optional Tukey HSD provides pairwise follow-up; clustering of significant rows describes patterns and is not an additional significance test. Scope: one factor (the condition) only; no covariates, batch terms, blocking or repeated measures, and no interaction tests.",
    "qc.cv": "CV = 100 * sample standard deviation / mean, within each condition, on linear intensities before imputation (log2 values are back-transformed first; when a normalization was applied, the normalized values are used, except after VSN). Sample SD uses ddof=1. Missing observations are excluded, not zero-filled; nonfinite CVs are excluded. This is descriptive variability, not a statistical significance test or a universal quality threshold.",
    "stats.volcano": "The x-axis is the selected contrast's log2 fold change; the first selected condition minus the second for the ordinary comparison, so a positive value means higher in condition 1. The control variant is (condition 1 - its control) minus (condition 2 - its control). The y-axis is -log10 of the selected p-metric, adjusted unless uncorrected values were selected. BH adjustment controls false discovery rate across the applicable test family. The selected test, pairing and thresholds determine inference; inspect results for their actual values. Log2FC=1 corresponds to a twofold ratio. Statistical significance does not prove biological importance or causality. Moderated test types (limma) shrink each row's variance towards a prior fitted to all rows and run in R; the robust variants fit that prior robustly. Fold changes and BH correction are unchanged; only the p-values differ. Scope: two-group comparisons and control-corrected contrasts only; pairing is the only blocking modelled, so there is no adjustment for batches, covariates or other factors.",
    "qc.coverage": "Counts observed entities in the selected workflow and samples. Peptide and precursor counts are different quantities. Detection depends on the declared scale: zero is non-detection on linear intensity but a valid value on log2 intensity. Missing values remain missing.",
    "qc.histogram": "Displays an intensity distribution. The selected aggregation and data source determine what each observation represents; do not infer individual-sample behavior from condition averages.",
    "qc.boxplot": "Summarizes intensity distributions with medians, quartiles and whiskers. A boxplot is descriptive and cannot establish batch effects or biological differences on its own.",
    "qc.pca": "Each dot represents one sample, positioned using its measurements across features. For PCA, the axes summarise variation between samples; the percentage on an axis is the share of total sample variation captured by that component, not how many percent the groups differ. PCA variance percentages do not apply to other projection methods. Scaling, source and method affect the picture. Separation alone does not establish significance or causality.",
    "qc.abundance": "Ranks abundance values. The selected condition and source determine the represented entities; a rank curve is not a differential-expression test.",
    "qc.correlation": "Pearson correlation matrix of all samples at once: every sample pair, so replicate agreement within and between conditions, the most and least similar pairs. Constant or insufficiently overlapping data can yield undefined correlations, which must not be invented or silently replaced with similarity.",
    "completeness.missingPlot": "Summarizes missing observations for the selected workflow and aggregation. Missingness does not establish biological absence or a missing-data mechanism.",
    "completeness.heatmap": "Displays missingness patterns. Aggregate missingness counts do not establish a particular row-level or block pattern without the corresponding matrix evidence.",
    "completeness.tables": "Reports completeness under the selected sample/condition rules. Totals and thresholds describe the chosen definition, not an automatic recommendation to filter or impute.",
    "stats.ora": (
        "Gene Ontology (GO) has three complementary branches: "
        "MF means Molecular Function, the molecular activity of a gene product (for example kinase activity, ATP binding or actin binding); "
        "BP means Biological Process, the broader process it contributes to (for example signal transduction, chromatin remodelling or an apoptotic process); "
        "CC means Cellular Component, its cellular location or complex (for example nucleus, membrane or cytoskeleton). "
        "A protein can have annotations in more than one branch, but is not guaranteed annotations in all three. "
        "These are vocabulary examples, not enriched findings in the current data. "
        "The app tests over-representation using a one-sided hypergeometric test and BH adjustment within each query's size-eligible, nonzero-overlap term family, "
        "across the eligible terms supplied from the selected database; the branch labels do not create separate correction families. "
        "Upregulated and downregulated query lists are analysed separately. The query and background must be appropriate. "
        "Display limits do not change adjusted p-values. Over-representation does not establish pathway activation, inhibition or causality; it describes how often annotations occur in the selected gene list relative to its background."
    ),
    "stats.pathwayHeatmap": "Displays selected pathway-associated entities on the selected scale. A heatmap pattern alone is not evidence of pathway activation or an enrichment test.",
    "singleProtein.boxplot": "Displays a selected protein or site's values by condition. Check identifier, source and any selected statistical test before interpreting a difference.",
    "singleProtein.lineplot": "Displays the selected entity across samples, retaining samples with missing measurements. Lines connect consecutive observed points. Show NA as low optionally draws missing positions at log2 0 (linear 1) without markers; these are display-only vertices, not measured or imputed values. Connecting points does not by itself establish a time course or paired design.",
    "singleProtein.heatmap": "Displays phosphorylation sites associated with the selected protein. Site-level abundance and protein-level abundance are distinct measurements.",
    "comparison.pearson": "Scatter plot of one selected sample or condition against another with its single Pearson r; not the all-sample matrix (that is QC Pipeline / Correlation Plot). Compares paired numerical observations using Pearson correlation. Correlation measures linear association, not agreement or causality; undefined coefficients require explicit reporting.",
    "comparison.venn": "Displays overlaps between selected entity sets. Set membership and identifier rules determine the counts; overlap alone is not an enrichment significance test.",
    "phospho.phosphositePlot": "Counts the phosphoproteins, phosphopeptides and phosphosites in the dataset after the localization cutoff (a bar per count; no single site is selected here, one site's values are on the Single Protein pages). Localization and abundance are different quantities; site regulation does not alone establish kinase activity.",
    "phospho.coverage": "Counts sites using the application's localization categories: Class I means localization probability at least 0.75; lower known probabilities are not Class I. Bar segments represent these site classes, not different conditions. Unknown localization is unclassified rather than assigned an invented quality class.",
    "phospho.distribution": "Summarizes the selected phosphosite distribution. Interpret the actual aggregation and localization selections from results rather than assuming defaults.",
    "phospho.sty": "Summarizes serine, threonine and tyrosine site composition. Composition reflects measured sites and experimental coverage, not unbiased proteome-wide phosphorylation prevalence.",
    "phospho.phosprotRegulation": "Tests every phosphosite with the selected volcano test and thresholds, then plots each protein with at least one significant site: the x-axis is the sum of its upregulated sites' log2 fold changes, the y-axis the absolute sum of its downregulated sites' log2 fold changes. It uses site-level data only; protein abundance is not included or corrected for. Inspect source and contrast direction; do not equate site abundance with occupancy.",
    "peptide.rtPlot": "Displays peptide retention-time information. Observed differences do not alone establish chromatography failure; check source, units and aggregation.",
    "peptide.modification": "Summarizes observed peptide modifications. Counts depend on peptide/precursor identity and selected quantity; they do not directly establish modification occupancy.",
    "peptide.missedCleavage": "A missed cleavage is an internal site where the selected protease could cut but the detected peptide still spans that site; it does not mean the whole protein remained undigested. Counts use distinct stripped peptide sequences per sample with finite, nonzero quantities. Trypsin cuts after lysine (K) or arginine (R), except before proline (P); Trypsin / Lys-C also permits the K-before-P cut. The expected terminal cut is not counted as missed. A descriptive distribution alone does not determine digestion quality without a stated criterion.",
    "summary.report": (
        "Summary / Report assembles the plots you select, with optional text above or below each, into one self-contained HTML file "
        "(Generate HTML Report, then Download HTML Report). The HTML report is a static document of the rendered results: it can be read "
        "and shared, but nothing in it can be rerun. Export notebook beside the download controls gives the reproducibility option: a ZIP "
        "with a Python notebook, requirements.txt, a manifest, the captured input tables and the settings, which recreates the report plots "
        "outside the app when its cells are run in order. Share the whole ZIP plus the original uploaded data file(s); set data_path in the notebook, use its required Python version and install requirements.txt. The captured inputs verify reproduction and do not replace the original upload. These files contain analysis data. The viewer is not an export from "
        "this page: it is a separate read-only build of CoPYlot (DataViewer launcher) that loads a prepared project from viewer_config.json "
        "so others can browse the same results in the app without changing them; it is set up by whoever packages the app (docs/INSTALLATION.md) "
        "and reruns nothing. To let someone rerun the analysis, give them the notebook ZIP and required original input files; to let them look, the HTML report or a viewer build."
    ),
    "summary.tables": "Summary / Tables lists the result tables of the session for download (Download CSV, which Excel opens). The tables are the ones the analysis pages produced with their current settings; nothing is recalculated here.",
    "summary.text": "Summary / Text holds the free-text entries used above and below the plots of the HTML report; it is editorial text, not a calculated result.",
    "peptide.sequenceCoverage": "Maps observed peptide sequences onto a protein sequence. Sequence coverage is not intensity coverage and does not prove absence of unobserved sequence regions.",
}


_CV_SOURCE_SENTENCES = {
    "none": "Use original uploaded linear intensities before application of imputation; back-transform log2 uploads first.",
    "vsn": "After VSN the unnormalized intensities of the same rows are used, because VSN's glog output back-transformed is not a linear intensity; imputation is not applied.",
}


def cv_method_text(normalization: str | None) -> str:
    """The CV method sentence for the normalization the session applied.

    With no normalization the original wording stands (uploaded linear
    intensities before imputation); a normalization is mentioned only when one
    was applied, VSN with its exception.
    """
    mode = str(normalization or "none").strip().lower()
    source = _CV_SOURCE_SENTENCES.get(mode) or (
        f"Computed on the {mode}-normalized values back-transformed to linear intensities, before imputation."
    )
    return (
        "CV = 100 * sample standard deviation / mean, within each condition. " + source
        + " Sample SD uses ddof=1. Missing observations are excluded, not zero-filled; nonfinite CVs are excluded."
        " This is descriptive variability, not a statistical significance test or a universal quality threshold."
    )


def go_database_sentence(organism: str = "human") -> str:
    """Which Gene Ontology release the bundled term table was built from, for a methods section.

    Read data/db/BUILD.log at call time, so a rebuilt table reports its own release.
    """
    from app.services.reference_db import reference_release

    gaf = {"human": "goa_human.gaf.gz", "mouse": "mgi.gaf.gz"}.get(str(organism).strip().lower(), "goa_human.gaf.gz")
    annotations = reference_release(gaf)
    ontology = reference_release("go-basic.obo.gz")
    if not annotations and not ontology:
        return ("The Gene Ontology release of the bundled term table is not recorded in this installation "
                "(see data/db/README.md); state the release date of the annotations you used when citing.")
    return (f"The bundled Gene Ontology table was built from GO release {annotations or 'unknown'} "
            f"(annotations {gaf.removesuffix('.gz')}; ontology go-basic.obo release {ontology or 'unknown'}), "
            "see data/db/README.md. Cite the Gene Ontology Consortium with that release date: "
            "Gene Ontology Consortium, Genetics 224, iyad031 (2023) and Ashburner et al., Nat Genet 25, 25-29 (2000); "
            "the current wording is at geneontology.org/docs/go-citation-policy.")


# The release is fixed for an installation, so the reviewed method text carries it
# wherever METHODS is read (answer catalogue, planner guidance, answer checker).
for _module in ("stats.ora", "stats.pathwayHeatmap"):
    METHODS[_module] = METHODS[_module] + " " + go_database_sentence()


def method_evidence(module: str) -> Evidence:
    if module == "workflow.overview":
        return Evidence(tool="method:workflow.overview", status="available", text=WORKFLOW_GUIDE)
    text = METHODS.get(module)
    return Evidence(tool=f"method:{module}", status="available" if text else "unavailable",
                    text=text or "No reviewed standalone method description exists for this module.")


def canvas_workflow_guidance() -> Evidence:
    return Evidence(tool="workflow:analysis.canvas", status="available", text=(
        "To build a focused canvas, ask for a protein or an exact phosphosite identifier—for example, "
        "‘Show everything on the canvas for phosphosite YOUR_SITE’. To expand a protein into sites, "
        "explicitly ask ‘Show all phosphosites on protein YOUR_PROTEIN on the canvas’. "
        "Protein, phosphosite and phosphoprotein data stay separate; the requested level must be loaded.\n\n"
        "You can specify the reference and target conditions. Otherwise the linked comparison is used, "
        "or the two available conditions; ambiguous choices require clarification. The app creates a "
        "highlighted volcano, an individual boxplot per feature and a shared lineplot. It keeps unrelated "
        "tiles and reuses the same focused group on repeated requests. To change an existing canvas's "
        "data level, save any work you need, then use an empty canvas.\n\n"
        "For manual setup, open Analysis, select a plot type and use + Add Plot. Each tile has its own "
        "controls and Plot styling section. Focus does not run enrichment or change preprocessing. "
        "Wait for plots to load before interpreting them; a feature excluded from the volcano comparison "
        "can still appear in its boxplot and lineplot. This is guidance only; nothing was changed."
    ))
