#!/usr/bin/env Rscript
# S7 Table: runs the four analysis steps that compute sheets but draw nothing,
# folds their CSVs and the ones F06.R and S7.R leave in c_data into one
# workbook, then deletes the CSVs. Run after F06.R and S7.R.

setwd(here::here())

source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F06/c_data"
needed <- file.path(DAT, c("module_grid/module_grid_summary.csv",
                           "panel_B_full_screen_bh.csv"))
if (!all(file.exists(needed))) {
  stop("no module grid or screen CSVs in ", DAT, "; run F06.R and S7.R first")
}

# Sourced into this session in their original order. They share its
# environment on purpose: _supp_loso_wgcna_refit.R swaps cor() for WGCNA's and
# restores it through the global environment.
# Writes classifier_pilot_*.csv
source("04_Figures/F06/a_script/_supp_prepare_roc.R")
# Writes panel_A_*.csv
source("04_Figures/F06/a_script/_supp_multivariate.R")
# Writes loso_auc/loso_auc_summary.csv
source("04_Figures/F06/a_script/_supp_loso_sensitivity.R")
# About 255 s; writes loso_auc/loso_wgcna_refit_*.csv
source("04_Figures/F06/a_script/_supp_loso_wgcna_refit.R")

# Restore the path _supp_multivariate.R points at F05's c_data.
DAT <- "04_Figures/F06/c_data"

f07 <- function(p) file.path(DAT, p)

message("F06 supplementary workbook")
f07_specs <- list(
  list(name = "module_grid_summary", path = f07("module_grid/module_grid_summary.csv")),
  list(name = "module_grid_curves", path = f07("module_grid/module_grid_curves.csv")),
  list(name = "panel_A_classifier_auc", path = f07("panel_A_multi_classifier_auc.csv")),
  list(name = "panel_A_feature_stability", path = f07("panel_A_feature_stability.csv")),
  list(name = "panel_A_permutation", path = f07("panel_A_permutation.csv")),
  list(name = "panel_A_roc_curves", path = f07("panel_A_roc_curves.csv")),
  list(name = "classifier_pilot_summary", path = f07("classifier_pilot_summary.csv")),
  list(name = "classifier_pilot_curves", path = f07("classifier_pilot_curves.csv")),
  list(name = "panel_B_full_screen", path = f07("panel_B_full_screen_bh.csv")),
  list(name = "loso_auc_summary", path = f07("loso_auc/loso_auc_summary.csv")),
  list(name = "loso_wgcna_refit_summary", path = f07("loso_auc/loso_wgcna_refit_summary.csv")),
  list(name = "loso_wgcna_refit_mod_stability", path = f07("loso_auc/loso_wgcna_refit_module_stability.csv"))
)
build_workbook(
  f07("F06_data.xlsx"),
  title = "S7 Table \u2014 module discrimination and phenotype coupling",
  description = "Phenotype-prediction outputs: univariate module-outcome ROCs, multivariate classifiers, LOSO cross-validation, age-stratified module-phenotype coupling, and per-module Jaccard stability.",
  overview_df = data.frame(
    Sheet = c(
      "module_grid_summary",
      "module_grid_curves",
      "panel_A_classifier_auc",
      "panel_A_feature_stability",
      "panel_A_permutation",
      "panel_A_roc_curves",
      "classifier_pilot_summary",
      "classifier_pilot_curves",
      "panel_B_full_screen",
      "loso_auc_summary",
      "loso_wgcna_refit_summary",
      "loso_wgcna_refit_mod_stability"
    ),
    Description = c(
      "Per-module univariate ROC summary: AUC, permutation p, BH q for each module-outcome pair",
      "Per-module ROC curves: TPR/FPR coordinates for plotting",
      "Multivariate classifier AUC: raw-protein (k=10), phenotype-only, ME-stack, and delta-ME classifiers",
      "Feature stability: top-10 raw-protein classifier features across LOSO folds",
      "Permutation null: observed vs permuted AUC distributions for each classifier",
      "Multivariate ROC curves: TPR/FPR coordinates for each classifier",
      "Classifier pilot summary: module eigengene ROCs for age discrimination (top modules)",
      "Classifier pilot curves: TPR/FPR coordinates for classifier pilot ROCs",
      "Panel B full screen: 180-test module-phenotype correlations (9 modules \u00d7 2 sources \u00d7 5 outcomes \u00d7 2 strata) with BH correction",
      "LOSO fixed-module: leave-one-subject-out AUCs using fixed full-sample module definitions",
      "LOSO + WGCNA refit: full-refit LOSO \u2014 top-12 pair AUCs with network refit on n\u22121 subjects per fold, training modules matched to full-sample by Jaccard",
      "LOSO module stability: per-full-sample-module mean/min Jaccard of training-fold vs full-sample assignments, + count of folds where best Jaccard fell below 0.5"
    ),
    stringsAsFactors = FALSE
  ),
  sheet_specs = f07_specs
)
cleanup_after_workbook(f07_specs,
  extra_subdirs = c(f07("loso_auc"), f07("module_grid"))
)

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("F06 complete")
