# S4b heatmap for one significance criterion, shared by panels A-D, on
# the shared functional-classification engine with
# cluster_within_group = TRUE: rows inside each group follow hierarchical
# clustering (Euclidean/Ward.D2 on row z-scores), not log2FC, because these
# panels show how BH-FDR and the Xiao pi-score diverge on protein selection.
# Each panel's off-criterion column (Pi status on FDR-gated panels, FDR status
# on Pi-gated ones) shows the overlap directly. min_group_size is smaller than
# the engine was first tuned for because these protein sets are smaller.

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot,
  ComplexHeatmap, circlize, cluster
)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

RPT <- "04_Figures/F03/b_reports/supp/panels"
DAT_HT <- "04_Figures/F03/c_data/sig_heatmaps"
for (d in c(RPT, DAT_HT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

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

# The letters are the caption's. S4b ships as four pages and the legend calls
# them (A) to (D); without them on the plot the reader has only page order.
build_sig_heatmap <- function(protein_df, logfc_cols, pi_col, fdr_col,
                              min_group_size, title, subtitle_fmt,
                              file_stub, name, comp_h = NULL) {
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
    rpt_png = RPT, rpt_pdf = RPT, dat_ht = DAT_HT,
    comp_w = 178, comp_h = comp_h,
    row_name_cap = Inf, row_name_pt = SIG_LABEL_PT
  )
  source("04_Figures/shared/comparison_panels/panel_heatmap_classified.R", local = TRUE)

  for (ext in c("png", "pdf")) {
    file.rename(
      file.path(RPT, sprintf("MAIN_panel_%s.%s", file_stub, ext)),
      file.path(RPT, sprintf("%s.%s", name, ext))
    )
  }

  # The composite saves each heatmap as one page of S4b at this height.
  attr(pH_heat, "height_mm") <- comp_h
  pH_heat
}
