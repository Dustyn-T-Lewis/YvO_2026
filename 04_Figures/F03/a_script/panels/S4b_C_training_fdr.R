#!/usr/bin/env Rscript
# S4b Figure C: clustered heatmap of training in younger adults, FDR < 0.05.

setwd(here::here())
source("04_Figures/F03/a_script/panels/_sig_heatmap.R", local = TRUE)

training_fdr <- dep |> filter(!is.na(adj.P.Val_Training_Young), adj.P.Val_Training_Young < 0.05)
invisible(build_sig_heatmap(
  protein_df = training_fdr,
  logfc_cols = c(`log2FC Tr. (Y)` = "logFC_Training_Young", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = "pi_score_Training_Young", fdr_col = NULL,
  min_group_size = 3,
  title = "(C) Training (Young) signature: FDR-significant (functional classification)",
  subtitle_fmt = "%d Training(Y) FDR < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "training_fdr_heatmap",
  name = "S4b_C_training_fdr"
))
