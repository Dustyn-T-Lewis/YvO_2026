#!/usr/bin/env Rscript
# S1a Figure B: per-protein missingness.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

prot_miss <- rowSums(is.na(int_norm$data_pre_outlier)) /
  ncol(int_norm$data_pre_outlier) * 100
miss_hist_df <- tibble(pct_miss = prot_miss)
q99 <- quantile(prot_miss, 0.99, na.rm = TRUE)

pB <- ggplot(miss_hist_df, aes(pct_miss)) +
  geom_histogram(binwidth = 2.5, fill = "#4393C3", color = "white", alpha = 0.85) +
  geom_vline(xintercept = q99, linetype = "dashed", color = "#B2182B", linewidth = 0.4) +
  annotate("text",
    x = q99, y = Inf, label = sprintf("99th = %.1f%%", q99),
    vjust = 1.6, hjust = 1.05, color = "#B2182B", size = 2.2
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    x = "% missing", y = "Proteins", tag = "B",
    title = "Per-protein missingness distribution",
    subtitle = sprintf(
      "%s proteins | median = %.1f%%",
      comma(nrow(miss_hist_df)), median(prot_miss)
    )
  ) +
  FIG_THEME

save_panel(pB, "S1a_B_protein_missingness", "panel_B", miss_hist_df)

invisible(pB)
