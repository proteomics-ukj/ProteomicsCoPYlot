# Proteomics CoPYlot: HEK293, WT versus GRK
functions_path <- Sys.getenv("PROTEOMICS_COPYLOT_FUNCTIONS", unset = "")
if (!nzchar(functions_path)) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  start <- if (!is.na(script)) dirname(normalizePath(script)) else getwd()
  candidates <- file.path(c(start, dirname(start), getwd(), dirname(getwd()), dirname(dirname(getwd()))),
                          "proteomicscopylot", "notebook_functions.R")
  functions_path <- candidates[file.exists(candidates)][1]
}
if (is.na(functions_path) || !file.exists(functions_path)) {
  stop("Could not find notebook_functions.R; set PROTEOMICS_COPYLOT_FUNCTIONS.", call. = FALSE)
}
source(functions_path)

# Load the example ---------------------------------------------------------
data <- load_example_data("protein")
meta <- load_example_data("protein_metadata")
sites <- load_example_data("phospho")
meta_sites <- load_example_data("phospho_metadata")
conditions <- c("WT", "GRK")
cat(nrow(data), "protein rows; replicates per condition:\n")
print(table(meta$condition))

# Filter proteins ----------------------------------------------------------
filtered <- filter_data(data, meta, 3, mode = "per_group")
filter_protein_contaminants <- TRUE
contaminant_categories <- 1:5
filtered <- filter_contaminants(filtered, enabled = filter_protein_contaminants,
                           categories = contaminant_categories, match_accessions = TRUE)
if (!nrow(filtered)) stop("No proteins remain; revise the filters.")
cat(nrow(data), "proteins,", nrow(filtered), "after filtering\n")

# Protein: completeness ----------------------------------------------------
fig <- missing_value_plot(filtered, meta, data_scale = "log2")
missing_features <- missing_protein_table(filtered, meta, id=FALSE)
print(head(missing_features, 10))
fig <- missing_value_heatmap(filtered, meta, id = FALSE, data_scale = "log2")
print(missing_sample_table(filtered, meta, id = FALSE, data_scale = "log2"))

# Imputation ---------------------------------------------------------------
imputed <- impute_values(filtered, meta, seed = 1)
fig <- impute_values(filtered, meta, ret = 1, seed = 1)
fig <- impute_values(filtered, meta, ret = 2, seed = 1)
fig <- impute_values(filtered, meta, ret = 3, seed = 1)
print(head(imputed, 10))

# QC -------------------------------------------------------------------------
linear <- inverse_log2_transform_data(filtered, meta)
fig <- coverage_plot(filtered, meta, id = FALSE, data_scale = "log2")
coverage <- qc_coverage_table(filtered, meta)
print(coverage)
fig <- boxplot_int(filtered, meta, data_scale = "log2")
intensity_summary <- qc_boxplot_table(filtered, meta)
print(head(intensity_summary, 10))
fig <- histo_int(filtered, meta, data_scale = "log2")
fig <- cov_plot(linear, meta, data_scale = "linear")
cv_summary <- qc_cv_table(linear, meta)
print(head(cv_summary, 10))
feature_cvs <- qc_feature_cvs(linear, meta)
print(head(feature_cvs, 10))
fig <- abundance_plot(linear, meta, workflow = "Protein", data_scale = "linear")
fig <- corr_plot(filtered, meta, id = FALSE, data_scale = "log2")
correlations <- sample_correlation_matrix(filtered, meta, id=FALSE)
print(head(correlations, 10))
fig <- pca_plot(filtered, meta)
scores <- pca_scores(filtered, meta)
print(head(scores, 10))

# Differential abundance -----------------------------------------------------
volcano <- volcano_data_f(filtered, meta, conditions[1], conditions[2], workflow = "Protein", test_type = "unpaired")
print(table(volcano$significance))
fig <- volcano_plot(volcano, conditions[1], conditions[2])
print(head(volcano, 10))
protein_untested <- volcano_untested_table(filtered, meta, conditions[1], conditions[2], kind="protein")
print(head(protein_untested, 10))
welch <- volcano_data_f(filtered, meta, conditions[1], conditions[2], workflow = "Protein", test_type = "unpaired_welch")
print(table(welch$significance))
volcano_imputed <- volcano_data_f(imputed, meta, conditions[1], conditions[2], workflow = "Protein", observed = filtered)
print(head(volcano_imputed[c("ProteinNames", "log2FC", "adj_pval", "imputed_values")], 5))

# ANOVA ------------------------------------------------------------------------
anova <- anova_data_f(filtered, meta)
fig <- withCallingHandlers(anova_plot(anova, meta), warning = function(w) {
  if (grepl("features dropped for clustering", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
})
print(head(anova, 10))

# GO over-representation -----------------------------------------------------
genes <- different_genes(filtered, meta, conditions[1], conditions[2], workflow = "Protein")
significant <- c(genes$Upregulated, genes$Downregulated)
ora <- enrichment_data(significant, background = genes$Background)
fig <- enrichment_analysis(significant, background = genes$Background)
print(head(ora, 10))

# Single protein and comparison --------------------------------------------
protein_id <- as.character(na.omit(filtered$ProteinNames)[1])
by_protein <- filtered
rownames(by_protein) <- by_protein$ProteinNames
fig <- boxplot_int_single_prot(by_protein, meta, protein_id, dots = TRUE)
protein_box <- single_protein_boxplot_table(filtered, meta, "protein", protein_id, unique(meta$condition))
print(head(protein_box, 10))
fig <- compare_prot_line(filtered, meta, conditions, c(protein_id), workflow = "Protein", id = FALSE)
protein_line <- single_protein_lineplot_table(filtered, meta, "protein", protein_id, conditions)
print(head(protein_line, 10))
fig <- correlation_xy_intensity(filtered, meta, conditions[1], conditions[2])
pearson <- pearson_table(filtered, meta, "protein", mode="Condition", condition1=conditions[1], condition2=conditions[2])
print(head(pearson, 10))
fig <- overlap_plot_condition(filtered, meta, conditions[1], conditions[2], data_scale = "log2")
overlap <- venn_table(filtered, meta, "protein", mode="Condition", first=conditions[1], second=conditions[2])
print(head(overlap, 10))

# Phosphosites ---------------------------------------------------------------
cat(nrow(sites), "phosphosite rows\n")
sites_filtered <- filter_data(sites, meta_sites, 3, mode = "per_group")
filter_phospho_contaminants <- TRUE
sites_filtered <- filter_contaminants(sites_filtered, enabled = filter_phospho_contaminants,
                                     categories = 1:5)
sites_filtered <- filter_localization(sites_filtered, 0)
fig <- simple_phos_site_plot(sites_filtered)
site_counts <- phospho_counts_table(sites_filtered)
print(head(site_counts, 10))
fig <- phossite_coverage_plot(sites_filtered, meta_sites, id = FALSE)
site_coverage <- phospho_coverage_table(sites_filtered, meta_sites, is_log2=TRUE)
print(head(site_coverage, 10))
fig <- sty_plot(sites_filtered)
site_residues <- phospho_sty_table(sites_filtered)
print(head(site_residues, 10))
fig <- missing_value_plot(sites_filtered, meta_sites, data_scale = "log2")
missing_sites <- missing_protein_table(sites_filtered, meta_sites, id=FALSE)
print(head(missing_sites, 10))
site_volcano <- volcano_data_f(sites_filtered, meta_sites, conditions[1], conditions[2], workflow = "Phosphosite")
fig <- volcano_plot(site_volcano, conditions[1], conditions[2])
print(head(site_volcano, 10))
phospho_untested <- volcano_untested_table(sites_filtered, meta_sites, conditions[1], conditions[2], kind="phospho")
print(head(phospho_untested, 10))
site_id <- as.character(na.omit(sites_filtered$PTM_Collapse_key)[1])
by_site <- sites_filtered
rownames(by_site) <- by_site$PTM_Collapse_key
fig <- boxplot_int_single_prot(by_site, meta_sites, site_id, dots = TRUE)
site_box <- single_protein_boxplot_table(sites_filtered, meta_sites, "phospho", site_id, unique(meta_sites$condition))
print(head(site_box, 10))
protein_of_site <- split_tokens(as.character(na.omit(sites_filtered$Protein_group)[1]))[1]
fig <- phossite_on_prot_regulation(sites_filtered, meta_sites, protein_of_site, conditions = conditions, id = FALSE)
site_heatmap <- single_protein_heatmap_table(sites_filtered, meta_sites, "phospho", protein_of_site, conditions, filter_m1=FALSE)
print(head(site_heatmap, 10))

# Phosphoproteins ------------------------------------------------------------
phosprot <- transform_phosprot(sites_filtered, meta_sites, mode = "robust_site_median",
                              data_scale = "log2", coverage_data = sites_filtered)
phosprot_filtered <- filter_data(phosprot, meta_sites, 3, mode = "per_group")
phosprot_filtered <- filter_contaminants(phosprot_filtered, enabled = filter_phospho_contaminants,
                                        categories = 1:5)
cat(nrow(sites_filtered), "retained sites;", nrow(phosprot_filtered), "phosphoproteins after filtering\n")
fig <- phossite_dist(phosprot_filtered, cutoff = 0.01)
site_distribution <- phospho_distribution_table(phosprot_filtered, cutoff=0.01)
print(site_distribution)
phosprot_volcano <- volcano_data_f(phosprot_filtered, meta_sites, conditions[1], conditions[2], workflow = "Phosphoprotein")
fig <- volcano_plot(phosprot_volcano, conditions[1], conditions[2])
print(head(phosprot_volcano, 10))
phosprot_untested <- volcano_untested_table(phosprot_filtered, meta_sites, conditions[1], conditions[2], kind="phosprot")
print(head(phosprot_untested, 10))
phosprot_id <- as.character(na.omit(phosprot_filtered$Phosphoprotein)[1])
fig <- compare_prot_line(phosprot_filtered, meta_sites, conditions, c(phosprot_id), workflow = "Phosphoprotein", id = FALSE)
phosprot_line <- single_protein_lineplot_table(phosprot_filtered, meta_sites, "phosprot", phosprot_id, conditions)
print(head(phosprot_line, 10))
