#!/usr/bin/env Rscript
# F03 Supp — four significance-criterion heatmaps: Aging x {FDR, Pi} and
# Training(Young) x {FDR, Pi}. Same functional-classification engine as F04
# panel D / F05 panel C, but with cluster_within_group = TRUE so rows inside
# each functional group are ordered by real hierarchical clustering
# (Euclidean/Ward.D2 on the row-z-scored matrix), not by log2FC -- these
# panels exist to show how the two significance criteria (BH-FDR vs. the
# Xiao pi-score) diverge on protein selection, not to re-litigate F04/F05's
# grouping story. Each panel's off-criterion column (FDR-gated panels show
# Pi status, Pi-gated panels show FDR status) makes the overlap between the
# two criteria visible directly on the plot. Mirrors 03_Figures/F03 on
# limma's own engine (03_DEP data) -- min_group_size values differ from
# smaller protein sets than the design was first drawn against.

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot,
  ComplexHeatmap, circlize, cluster
)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

BASE <- "04_Figures/F03"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT_HT <- file.path(BASE, "c_data", "sig_heatmaps")
for (d in c(PNL_PNG, PNL_PDF, DAT_HT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
id_cols <- c("uniprot_id", "protein", "gene", "description")
sample_cols <- setdiff(names(dep), id_cols)
sample_cols <- sample_cols[grepl("_(Pre|Post)$", sample_cols)]
universe <- unique(dep$gene[!is.na(dep$gene)])

# Every protein is named on these panels, so the page grows with the list
# rather than the list shrinking to the page. 1.3 line-heights per row at the
# label size, plus the header, legend and column labels above and below.
SIG_LABEL_PT <- 5
sig_heatmap_height <- function(n_rows, label_pt = SIG_LABEL_PT, furniture_mm = 58) {
  furniture_mm + n_rows * label_pt * 1.3 * 25.4 / 72
}

build_sig_heatmap <- function(protein_df, logfc_cols, pi_col, fdr_col,
                              min_group_size, title, subtitle_fmt,
                              file_stub, comp_h = NULL) {
  if (is.null(comp_h)) comp_h <- sig_heatmap_height(nrow(protein_df))
  cfg <- list(
    fig_id = "F03",
    protein_df = protein_df,
    sample_cols = sample_cols,
    logfc_cols = logfc_cols,
    pi_col = pi_col,
    fdr_col = fdr_col,
    min_group_size = min_group_size,
    cluster_within_group = TRUE,
    universe_genes = universe,
    title = title,
    subtitle_fmt = subtitle_fmt,
    file_stub = file_stub,
    rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat_ht = DAT_HT,
    comp_w = 178, comp_h = comp_h,
    row_name_cap = Inf, row_name_pt = SIG_LABEL_PT
  )
  source("04_Figures/shared/comparison_panels/panel_heatmap_classified.R", local = TRUE)

  file.rename(
    file.path(PNL_PNG, sprintf("MAIN_panel_%s.png", file_stub)),
    file.path(PNL_PNG, sprintf("SUPP_panel_%s.png", file_stub))
  )
  file.rename(
    file.path(PNL_PDF, sprintf("MAIN_panel_%s.pdf", file_stub)),
    file.path(PNL_PDF, sprintf("SUPP_panel_%s.pdf", file_stub))
  )

  ggsave(file.path(RPT_PNG, sprintf("SUPP_F03_%s.png", file_stub)), pH_heat,
    width = 178, height = comp_h, units = "mm", dpi = 300, bg = "white",
    device = ragg::agg_png
  )
  ggsave(file.path(RPT_PDF, sprintf("SUPP_F03_%s.pdf", file_stub)), pH_heat,
    width = 178, height = comp_h, units = "mm",
    device = get_pdf_device(), bg = "white"
  )
}

# The letters are the caption's. S4b ships as four pages and the legend calls
# them (A) to (D); without them on the plot the reader has only page order.
aging_fdr <- dep |> filter(!is.na(adj.P.Val_Aging), adj.P.Val_Aging < 0.05)
build_sig_heatmap(
  protein_df = aging_fdr,
  logfc_cols = c(`log2FC Aging` = "logFC_Aging", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = "pi_score_Aging", fdr_col = NULL,
  min_group_size = 8,
  title = "(A) Aging signature: FDR-significant (functional classification)",
  subtitle_fmt = "%d Aging FDR < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "aging_fdr_heatmap"
)

aging_pi <- dep |> filter(!is.na(pi_score_Aging), pi_score_Aging < 0.05)
build_sig_heatmap(
  protein_df = aging_pi,
  logfc_cols = c(`log2FC Aging` = "logFC_Aging", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = NULL, fdr_col = "adj.P.Val_Aging",
  min_group_size = 6,
  title = "(B) Aging signature: Pi-significant (functional classification)",
  subtitle_fmt = "%d Aging Π < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "aging_pi_heatmap"
)

training_fdr <- dep |> filter(!is.na(adj.P.Val_Training_Young), adj.P.Val_Training_Young < 0.05)
build_sig_heatmap(
  protein_df = training_fdr,
  logfc_cols = c(`log2FC Tr. (Y)` = "logFC_Training_Young", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = "pi_score_Training_Young", fdr_col = NULL,
  min_group_size = 3,
  title = "(C) Training (Young) signature: FDR-significant (functional classification)",
  subtitle_fmt = "%d Training(Y) FDR < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "training_fdr_heatmap"
)

training_pi <- dep |> filter(!is.na(pi_score_Training_Young), pi_score_Training_Young < 0.05)
build_sig_heatmap(
  protein_df = training_pi,
  logfc_cols = c(`log2FC Tr. (Y)` = "logFC_Training_Young", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = NULL, fdr_col = "adj.P.Val_Training_Young",
  min_group_size = 3,
  title = "(D) Training (Young) signature: Pi-significant (functional classification)",
  subtitle_fmt = "%d Training(Y) Π < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "training_pi_heatmap"
)

message("F03 supp: 4 significance-criterion heatmaps done (Aging FDR/Pi, Training(Y) FDR/Pi)")
