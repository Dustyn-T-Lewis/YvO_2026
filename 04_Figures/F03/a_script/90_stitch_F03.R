#!/usr/bin/env Rscript
# F03 — Volcano Rings: Master Orchestrator

setwd(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F03/c_data"

source("04_Figures/F03/a_script/02_supp_panels.R")
source("04_Figures/F03/a_script/01_main_panels.R")
source("04_Figures/F03/a_script/_supp_sig_heatmaps.R")

# Build xlsx from ring terms + significance-heatmap groupings + supp CSVs.
#
# The four per-contrast limma tables used to be copied in here as MAIN_<ctr>.
# They are S10 Table's Aging, Training_Young, Training_Old and Interaction
# sheets byte for byte, so the package was shipping 8,424 rows x 15 columns
# twice and a reader had no way to tell which copy was authoritative. S10 Table
# is the record; this workbook carries what the figure draws.

f03_specs <- list()
for (tag in c("A", "B", "C", "D")) {
  ring_path <- file.path(DAT, paste0("panel_", tag), "ring_terms.csv")
  if (file.exists(ring_path)) {
    f03_specs <- c(f03_specs, list(list(name = paste0("RING_", tag), path = ring_path)))
  }
}
# Excel sheet names cap at 31 chars; "SIGHEAT_<file_stub>_cluster_classification"
# blows past that, so short codes stand in for the file_stub.
sig_stub_codes <- c(
  aging_fdr_heatmap = "AgeFDR", aging_pi_heatmap = "AgePi",
  training_fdr_heatmap = "TrFDR", training_pi_heatmap = "TrPi"
)
for (p in list.files(file.path(DAT, "sig_heatmaps"), pattern = "\\.csv$", full.names = TRUE)) {
  base <- tools::file_path_sans_ext(basename(p))
  suffix <- if (grepl("_cluster_classification$", base)) "_cls" else "_grp"
  stub <- sub("_(groups|cluster_classification)$", "", base)
  f03_specs <- c(f03_specs, list(list(name = paste0(sig_stub_codes[[stub]], suffix), path = p)))
}
supp_csvs <- list.files(file.path(DAT, "supp"), pattern = "\\.csv$", full.names = TRUE)
for (p in supp_csvs) {
  f03_specs <- c(f03_specs, list(list(
    name = paste0("SUPP_", tools::file_path_sans_ext(basename(p))),
    path = p
  )))
}

# Keyed by sheet name, not built positionally: the heatmap and supp sheets come
# from globs, so their order follows the filesystem and an unlisted sheet has to
# fail loudly rather than inherit its neighbour's description.
f03_descriptions <- c(
  RING_A = "Panel A: the pathway terms the aging volcano ring labels",
  RING_B = "Panel B: pathway terms the training (younger) ring labels",
  RING_C = "Panel C: pathway terms the training (older) ring labels",
  RING_D = "Panel D: pathway terms the interaction ring labels",
  AgeFDR_grp = "S4b: group assignment, aging at FDR < 0.05 (278 proteins)",
  AgeFDR_cls = "S4b: per-protein cluster classification, aging at FDR < 0.05",
  AgePi_grp = "S4b: group assignment, aging at \u03a0 < 0.05 (195 proteins)",
  AgePi_cls = "S4b: per-protein cluster classification, aging at \u03a0 < 0.05",
  TrFDR_grp = "S4b: group assignment, training (younger) at FDR < 0.05 (135)",
  TrFDR_cls = "S4b: cluster classification, training (younger) at FDR < 0.05",
  TrPi_grp = "S4b: group assignment, training (younger) at \u03a0 < 0.05 (99)",
  TrPi_cls = "S4b: cluster classification, training (younger) at \u03a0 < 0.05",
  SUPP_panel_P_Value = "S4a panel A: raw p-values, all four contrasts, as drawn",
  SUPP_panel_pi_score = "S4a panel B: \u03a0 scores, all four contrasts",
  SUPP_panel_adj_P_Val = "S4a panel C: BH-adjusted p-values, all four contrasts",
  SUPP_panel_C_ma = "S4a panel D: MA-plot source, mean log2 intensity against log2 fold change"
)
f03_sheet_names <- vapply(f03_specs, `[[`, character(1), "name")
missing_desc <- setdiff(f03_sheet_names, names(f03_descriptions))
if (length(missing_desc)) {
  stop("no F03 description for: ", paste(missing_desc, collapse = ", "))
}

build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  title = "S4 Table \u2014 per-contrast results and significance criteria",
  description = "Source data for Figure 3, S4a Figure and S4b Figure: the pathway terms each volcano ring labels, the FDR and \u03a0 heatmap groupings, and the model diagnostics S4a draws. Per-contrast limma results are in S10 Table.",
  overview_df = data.frame(
    Sheet = f03_sheet_names,
    Description = unname(f03_descriptions[f03_sheet_names])
  ),
  sheet_specs = f03_specs
)
cleanup_after_workbook(f03_specs,
  extra_subdirs = c(
    file.path(DAT, "panel_A"), file.path(DAT, "panel_B"),
    file.path(DAT, "panel_C"), file.path(DAT, "panel_D"),
    file.path(DAT, "sig_heatmaps"), file.path(DAT, "supp")
  )
)

message("F03 complete")
