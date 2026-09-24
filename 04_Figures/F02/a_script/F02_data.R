#!/usr/bin/env Rscript
# S3 Table: folds the CSVs that F02.R and S3.R leave in c_data into one
# workbook, then deletes them. Run after those two.

setwd(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F02/c_data"

audit_csvs <- list.files(DAT, pattern = "^(audit_|panel_|SUPP_panel_)", full.names = TRUE)
if (!length(audit_csvs)) {
  stop("no panel CSVs in ", DAT, "; run F02.R and S3.R first")
}
specs <- lapply(audit_csvs, \(p) list(name = tools::file_path_sans_ext(basename(p)),
                                       path = p))

# Keyed by sheet name rather than positionally: the specs come from a glob, so
# their order follows the filesystem and an unlisted sheet has to fail loudly
# rather than inherit its neighbour's description.
f02_descriptions <- c(
  panel_A_pca_variance_ci = "Panel A: PCA variance explained per component, with bootstrap CI",
  panel_B_logfc_long = "Panel B: per-protein log2 fold change, all four contrasts in long form",
  panel_B_stats = "Panel B: distribution summary per contrast (median, IQR, tail counts)",
  panel_D_upset_overlap_enrich = "Panel D: UpSet intersection sizes, with ORA per intersection",
  panel_E_fgsea_counts = "Panel E: significant gene-set counts per contrast and collection",
  panel_E_fgsea_sig = "Panel E: the significant gene sets, with NES and adjusted p",
  panel_F_barcode_enrichment = "Panel F: barcode enrichment statistics per contrast",
  SUPP_panel_A_cv_delta = "S3a: change in per-protein CV% from Pre to Post, per age group",
  SUPP_panel_A_cv_scatter = "S3a: per-protein CV% Pre against Post, as the scatter draws it",
  SUPP_panel_B_median_cv_ci = "S3b: median inter-individual CV% per group, bootstrap 95% CI",
  SUPP_panel_B_wilcoxon = "S3b: Wilcoxon tests behind the CV comparisons in panel B",
  SUPP_panel_C_imputed = "S3c: per-subject intra-individual variability and imputed fraction",
  SUPP_panel_C_wilcoxon = "S3c: Wilcoxon tests behind the intra-individual comparisons"
)
sheet_names <- vapply(specs, `[[`, character(1), "name")
missing_desc <- setdiff(sheet_names, names(f02_descriptions))
if (length(missing_desc)) {
  stop("no F02 description for: ", paste(missing_desc, collapse = ", "))
}
build_workbook(
  file.path(DAT, "F02_data.xlsx"),
  title = "S3 Table \u2014 proteome overview and differential expression",
  description = "Source data for Figure 2 and S3 Figure: principal components, effect-size distributions, contrast overlaps, pathway counts and coefficient-of-variation diagnostics.",
  overview_df = data.frame(Sheet = sheet_names,
                           Description = unname(f02_descriptions[sheet_names])),
  sheet_specs = specs)
cleanup_after_workbook(specs)

message("F02 complete")
