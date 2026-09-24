#!/usr/bin/env Rscript
# S4b Figure D: clustered heatmap of training in younger adults, Pi < 0.05.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_sig_heatmap.R", local = TRUE)

training_pi <- dep |> filter(!is.na(pi_score_Training_Young), pi_score_Training_Young < 0.05)
invisible(build_sig_heatmap(
  protein_df = training_pi,
  logfc_cols = c(`log2FC Tr. (Y)` = "logFC_Training_Young", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = NULL, fdr_col = "adj.P.Val_Training_Young",
  min_group_size = 3,
  title = "(D) Training (Young) signature: Pi-significant (functional classification)",
  subtitle_fmt = "%d Training(Y) Π < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "training_pi_heatmap",
  name = "S4b_D_training_pi"
))
