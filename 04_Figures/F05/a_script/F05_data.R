#!/usr/bin/env Rscript
# S6 Table: builds c_data/F05_data.xlsx from the network outputs of
# YvO_WGCNA_run.R, the CSVs F05.R, S6.R and S8.R leave in c_data, the hub edge
# table and the preservation Z-summaries, then deletes the intermediates. Run
# after those three.

setwd(here::here())

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/figure_supplement_helpers.R")

pacman::p_load(dplyr, tidyr)

BASE <- "04_Figures/F05"

panel_csvs <- file.path(BASE, "c_data", c(
  "01_panel_A_heatmap_data.csv", "panel_B_module_fgsea.csv",
  "03_panel_B_heatmap_zscores.csv", "supp/a05_sft_fit_indices.csv",
  "asupp_B_QC_dendrogram_SUPP_data.csv", "supp/a03_compartment_enrichment.csv",
  "supp/a02_bicor_sensitivity.csv"
))
if (!all(file.exists(panel_csvs))) {
  stop(
    "missing panel CSVs in ", file.path(BASE, "c_data"),
    "; run F05.R, S6.R and S8.R first"
  )
}

# Fetches from STRING on first run and reads the cache afterwards, so the
# module labels' "densest STRING cluster" claim has something behind it that a
# reader can check. No network needed once the cache is in place.
source("04_Figures/shared/build_string_cluster_cache.R")

DAT <- file.path(BASE, "c_data")
f06 <- function(p) file.path(DAT, p)

# RDS matrix files converted to data frames for Excel embedding
.shared <- readRDS(f06("shared_objects.rds"))
.MEs <- matrix_to_df(as.matrix(readRDS(f06("MEs.rds"))), "sample_id")
.me_pre <- matrix_to_df(as.matrix(readRDS(f06("me_pre.rds"))), "subject_key")
.me_post <- matrix_to_df(as.matrix(readRDS(f06("me_post.rds"))), "subject_key")
.delta_me <- matrix_to_df(as.matrix(readRDS(f06("delta_me.rds"))), "subject_key")

# Consolidate redundant cohort/metric/mode/check splits into long sheets
# Folds 22 sheets into 5 long-format sheets, while preserving sheet names that
# downstream consumers (notably F06 panels and supp scripts) read directly.
.read_long <- function(path, cohort, metric) {
  df <- safe_read(path)
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  long <- tidyr::pivot_longer(df, -module, names_to = "trait", values_to = "value")
  long$cohort <- cohort
  long$metric <- metric
  long[, c("module", "trait", "cohort", "metric", "value")]
}
.collect <- function(stage) {
  variants <- list(
    list(c = "combined", m = "cor", p = sprintf("wgcna/wgcna_%s_trait_correlations.csv", stage)),
    list(c = "young", m = "cor", p = sprintf("wgcna/wgcna_%s_trait_correlations_young.csv", stage)),
    list(c = "old", m = "cor", p = sprintf("wgcna/wgcna_%s_trait_correlations_old.csv", stage)),
    list(c = "combined", m = "pval_bh", p = sprintf("wgcna/wgcna_%s_trait_pvalues_bh.csv", stage)),
    list(c = "young", m = "pval_bh", p = sprintf("wgcna/wgcna_%s_trait_pvalues_bh_young.csv", stage)),
    list(c = "old", m = "pval_bh", p = sprintf("wgcna/wgcna_%s_trait_pvalues_bh_old.csv", stage)),
    list(c = "young", m = "pval_raw", p = sprintf("wgcna/wgcna_%s_trait_pvalues_raw_young.csv", stage)),
    list(c = "old", m = "pval_raw", p = sprintf("wgcna/wgcna_%s_trait_pvalues_raw_old.csv", stage))
  )
  do.call(rbind, lapply(variants, function(v) .read_long(f06(v$p), v$c, v$m)))
}
.baseline_trait_assoc <- .collect("baseline")
.change_trait_assoc <- .collect("change")

# Module-trait associations (main combined design): combine cor + pval_bh
.mt_cor <- safe_read(f06("wgcna/wgcna_module_trait_correlations.csv"))
.mt_bh <- safe_read(f06("wgcna/wgcna_module_trait_pvalues_bh.csv"))
.module_trait_assoc <- rbind(
  cbind(tidyr::pivot_longer(.mt_cor, -module, names_to = "trait", values_to = "value"),
    metric = "cor"
  ),
  cbind(tidyr::pivot_longer(.mt_bh, -module, names_to = "trait", values_to = "value"),
    metric = "pval_bh"
  )
)
.module_trait_assoc <- .module_trait_assoc[, c("module", "trait", "metric", "value")]

# Module enrichment: strict + relaxed (+ mode)
.enr_strict <- safe_read(f06("wgcna/wgcna_module_enrichment.csv"))
.enr_relaxed <- safe_read(f06("03_panel_B_triptych_enrichment.csv"))
if (!is.null(.enr_strict)) .enr_strict$mode <- "strict"
if (!is.null(.enr_relaxed)) .enr_relaxed$mode <- "relaxed"
.module_enrichment <- dplyr::bind_rows(.enr_strict, .enr_relaxed)

# LMM diagnostics: contrast + stratified (+ check)
.lmm_c <- safe_read(f06("wgcna/wgcna_lmm_contrast_check.csv"))
.lmm_s <- safe_read(f06("wgcna/wgcna_lmm_stratified_check.csv"))
if (!is.null(.lmm_c)) {
  .lmm_c$check <- "contrast"
  .lmm_c$age_group <- NA_character_
}
if (!is.null(.lmm_s)) {
  .lmm_s$check <- "stratified"
  .lmm_s$contrast <- NA_character_
}
.lmm_diagnostics <- dplyr::bind_rows(.lmm_c, .lmm_s)

# mod_bio_labels.csv carries each module's top five ORA terms as one
# semicolon-joined evidence string, which no reader can sort or filter. Split it
# back into rows and re-attach each term's statistics from the enrichment run.
.ora_terms <- local({
  labs <- safe_read(f06("mod_bio_labels.csv"))
  enr <- safe_read(f06("wgcna/wgcna_module_enrichment.csv"))
  if (is.null(labs) || is.null(enr)) {
    return(NULL)
  }
  # A term name can be shared by two collections; the evidence string was built
  # from the best-scoring copy, so keep that one and the join stays 1:1.
  kept <- enr |>
    dplyr::filter(dedup_status == "kept") |>
    dplyr::arrange(module, padj, pathway) |>
    dplyr::distinct(module, Description, .keep_all = TRUE) |>
    dplyr::select(
      module_color = module, term = Description, pathway, database,
      padj, pval, overlap, size, foldEnrichment
    )
  terms <- labs |>
    dplyr::select(module_color, module_id, bio_label, ora_label, evidence) |>
    tidyr::separate_rows(evidence, sep = "; ") |>
    dplyr::rename(term = evidence) |>
    dplyr::group_by(module_color) |>
    dplyr::mutate(term_rank = dplyr::row_number()) |>
    dplyr::ungroup()
  out <- terms |>
    dplyr::left_join(kept, by = c("module_color", "term")) |>
    dplyr::arrange(match(module_color, labs$module_color), term_rank)
  stopifnot(
    "an ORA evidence term has no match in wgcna_module_enrichment.csv" =
      !anyNA(out$padj),
    "splitting the evidence string changed the term count" =
      nrow(out) == nrow(terms)
  )
  as.data.frame(out[, c(
    "module_color", "module_id", "bio_label", "ora_label", "term_rank",
    "term", "pathway", "database", "padj", "pval", "overlap", "size",
    "foldEnrichment"
  )])
})

source("04_Figures/F05/a_script/_module_labels.R")

f06_specs <- list(
  list(name = "panel_A_heatmap", path = f06("01_panel_A_heatmap_data.csv")),
  list(name = "panel_B_module_fgsea", path = f06("panel_B_module_fgsea.csv")),
  list(name = "WGCNA_module_assignments", path = f06("wgcna/wgcna_module_assignments.csv")),
  list(name = "WGCNA_mod_bio_labels", path = f06("mod_bio_labels.csv")),
  list(name = "WGCNA_module_core", df = .module_core),
  list(name = "WGCNA_module_ora", df = .module_ora),
  list(name = "WGCNA_string_clusters", path = f06("wgcna_string_clusters.csv")),
  list(name = "WGCNA_ora_terms", df = .ora_terms),
  list(name = "WGCNA_hub_proteins", path = f06("wgcna/wgcna_hub_proteins.csv")),
  list(name = "WGCNA_kme_all", path = f06("wgcna_kme_all.csv")),
  list(name = "WGCNA_hub_network_edges", path = f06("04_panel_D_hub_network.csv")),
  list(name = "WGCNA_protein_zscores_by_group", path = f06("03_panel_B_heatmap_zscores.csv")),
  list(name = "WGCNA_module_enrichment", df = .module_enrichment),
  list(name = "WGCNA_module_trait_assoc", df = .module_trait_assoc),
  list(name = "WGCNA_baseline_trait_assoc", df = .baseline_trait_assoc),
  list(name = "WGCNA_change_trait_assoc", df = .change_trait_assoc),
  list(name = "WGCNA_sft_summary", path = f06("wgcna/wgcna_sft_summary.csv")),
  list(name = "WGCNA_lmm_diagnostics", df = .lmm_diagnostics),
  list(name = "WGCNA_module_preservation", path = f06("05_panel_E_preservation.csv")),
  list(name = "WGCNA_gs_phenotype_choices", path = f06("wgcna/gs_phenotype_choices.csv")),
  list(name = "metadata_subj_age", path = f06("subj_age.csv")),
  list(name = "metadata_pheno_wide", path = f06("pheno_wide.csv")),
  list(name = "MEs", df = .MEs),
  list(name = "me_pre", df = .me_pre),
  list(name = "me_post", df = .me_post),
  list(name = "delta_me", df = .delta_me),
  list(name = "common_subj", df = data.frame(subject_key = .shared$common_subj, stringsAsFactors = FALSE)),
  list(name = "SUPP_panel_A_sft_fit", path = f06("supp/a05_sft_fit_indices.csv")),
  list(name = "SUPP_panel_B_dendrogram", path = f06("asupp_B_QC_dendrogram_SUPP_data.csv")),
  list(name = "SUPP_panel_C_compartment", path = f06("supp/a03_compartment_enrichment.csv")),
  list(name = "SUPP_panel_D_bicor_sensitivity", path = f06("supp/a02_bicor_sensitivity.csv"))
)

f06_overview <- data.frame(
  Sheet = c(
    "panel_A_heatmap", "panel_B_module_fgsea",
    "WGCNA_module_assignments", "WGCNA_mod_bio_labels",
    "WGCNA_hub_proteins", "WGCNA_hub_network_edges",
    "WGCNA_protein_zscores_by_group", "WGCNA_module_enrichment",
    "WGCNA_module_trait_assoc", "WGCNA_baseline_trait_assoc",
    "WGCNA_change_trait_assoc", "WGCNA_sft_summary",
    "WGCNA_lmm_diagnostics", "WGCNA_module_preservation",
    "WGCNA_gs_phenotype_choices", "metadata_subj_age", "metadata_pheno_wide",
    "MEs", "me_pre", "me_post", "delta_me", "common_subj",
    "WGCNA_module_core", "WGCNA_module_ora",
    "WGCNA_string_clusters", "WGCNA_ora_terms", "WGCNA_kme_all",
    "SUPP_panel_A_sft_fit", "SUPP_panel_B_dendrogram",
    "SUPP_panel_C_compartment", "SUPP_panel_D_bicor_sensitivity"
  ),
  Description = c(
    "Figure 5A: the module-trait heatmap as drawn, one row per module and contrast, with the correlation and its BH p-value",
    "Figure 5B: normalised enrichment score per module for each contrast",
    "Which module every protein was assigned to, grey meaning unassigned",
    "The two-part display name of each module, and the evidence for each half",
    "The ten highest-connectivity proteins in each module",
    "The co-expression edges drawn between those hub proteins",
    "Group-mean z-score for every protein, in the four age-by-timepoint cells",
    "Over-representation results for each module, across all gene-set collections",
    "Module against trait, as linear mixed-model contrasts with BH correction",
    "Module against trait at baseline only, cross-sectionally",
    "Change in module eigengene against change in outcome, per age stratum",
    "Scale-free topology fit at the chosen soft-thresholding power",
    "Convergence and degrees-of-freedom diagnostics for the mixed models",
    "Module preservation Z-summary between the pre and post networks",
    "Which phenotypes entered the module-trait screen, and why",
    "Participant identifier, age and age group",
    "Every phenotype measure, one row per participant",
    "Module eigengene for every sample",
    "Module eigengene at baseline, one row per participant",
    "Module eigengene after training, one row per participant",
    "The change in module eigengene, post minus pre",
    "The participants with both timepoints, which the paired analyses use",
    "Largest protein family in each module's kME >= 0.6 core",
    paste0(
      "GO:BP term on each module's label; it matched the core family in ",
      "seven of the nine modules and diverged in turquoise and black"
    ),
    "STRING MCL clusters per module (human, inflation 1.8, score >= 400)",
    "Top over-representation terms per module, one row per term with padj",
    "Module membership (kME) of every assigned protein against its own module eigengene, ranked within module",
    "Scale-free topology fit across candidate soft-thresholding powers",
    "Protein dendrogram and module colour assignment",
    "Subcellular compartment enrichment per module",
    "Module stability under bicor against Pearson correlation"
  )
)

# Two of the specs above read CSVs these scripts write, so they run before the
# workbook is built. Neither draws anything.
source("04_Figures/F05/a_script/_supp_mod_hub.R")
source("04_Figures/F05/a_script/_supp_preservation.R")

build_workbook(
  file.path("04_Figures/F05", "c_data", "F05_data.xlsx"),
  title = "S6 Table \u2014 co-expression network",
  description = "Source data for Figure 5, S6 Figure and S8 Figure: module assignments, eigengenes, module\u2013trait associations, the evidence behind each module name, and network construction diagnostics.",
  overview_df = f06_overview,
  sheet_specs = f06_specs
)
cleanup_after_workbook(f06_specs,
  extra_subdirs = character(0),
  preserve_patterns = c(
    "^00_input/", "^01_normalization/", "^02_imputation/", "^03_DEP/",
    "^04_Figures/shared/",
    "^04_Figures/F05/c_data/wgcna/",
    "^04_Figures/F05/c_data/.*\\.rds$",
    "^04_Figures/F05/c_data/key_modules\\.txt$",
    "^04_Figures/F05/c_data/mod_bio_labels\\.csv$",
    "^04_Figures/F05/c_data/wgcna_string_clusters\\.csv$",
    "^04_Figures/F05/c_data/wgcna_kme_all\\.csv$",
    "^04_Figures/F05/c_data/meta\\.csv$",
    "^04_Figures/F05/c_data/imp_annotations\\.csv$",
    "^04_Figures/F05/c_data/subj_age\\.csv$",
    "^04_Figures/F05/c_data/pheno_wide\\.csv$"
  )
)

message("F05 complete")
