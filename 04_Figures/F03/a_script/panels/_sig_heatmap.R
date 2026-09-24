# S4b heatmap for one significance criterion, shared by panels A-D. Same
# functional-classification engine as F04 panel D / F05 panel C, but with cluster_within_group = TRUE so rows inside
# each functional group are ordered by real hierarchical clustering
# (Euclidean/Ward.D2 on the row-z-scored matrix), not by log2FC -- these
# panels exist to show how the two significance criteria (BH-FDR vs. the
# Xiao pi-score) diverge on protein selection, not to re-litigate F04/F05's
# grouping story. Each panel's off-criterion column (FDR-gated panels show
# Pi status, Pi-gated panels show FDR status) makes the overlap between the
# two criteria visible directly on the plot. Mirrors 03_Figures/F03 on
# limma's own engine (03_DEP data) -- min_group_size values differ from
# smaller protein sets than the design was first drawn against.

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot,
  ComplexHeatmap, circlize, cluster
)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

RPT <- "04_Figures/F03/b_reports/panels"
DAT_HT <- "04_Figures/F03/c_data/sig_heatmaps"
for (d in c(RPT, DAT_HT)) {
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
