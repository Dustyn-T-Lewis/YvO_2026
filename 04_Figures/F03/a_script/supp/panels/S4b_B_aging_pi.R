#!/usr/bin/env Rscript
# S4b Figure B: clustered heatmap of aging, Pi < 0.05.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_sig_heatmap.R", local = TRUE)

aging_pi <- dep |> filter(!is.na(pi_score_Aging), pi_score_Aging < 0.05)
invisible(build_sig_heatmap(
  protein_df = aging_pi,
  logfc_cols = c(`log2FC Aging` = "logFC_Aging", `log2FC Tr. (O)` = "logFC_Training_Old"),
  pi_col = NULL, fdr_col = "adj.P.Val_Aging",
  min_group_size = 6,
  title = "(B) Aging signature: Pi-significant (functional classification)",
  subtitle_fmt = "%d Aging Π < 0.05 proteins | abundance per sample (row z-score) | hierarchical clustering within %d functional groups (GO Slim + Hallmark)",
  file_stub = "aging_pi_heatmap",
  name = "S4b_B_aging_pi"
))
