#!/usr/bin/env Rscript
# S4b Figure A: clustered heatmap of aging, FDR < 0.05.

setwd(here::here())
source("04_Figures/F03/a_script/panels/_sig_heatmap.R", local = TRUE)

aging_fdr <- dep |> filter(!is.na(adj.P.Val_Aging), adj.P.Val_Aging < 0.05)
invisible(build_sig_heatmap(
  protein_df = aging_fdr,
  logfc_cols = c(`log2FC Aging` = "logFC_Aging", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = "pi_score_Aging", fdr_col = NULL,
  min_group_size = 8,
  title = "(A) Aging signature: FDR-significant (functional classification)",
  subtitle_fmt = "%d Aging FDR < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "aging_fdr_heatmap",
  name = "S4b_A_aging_fdr"
))
