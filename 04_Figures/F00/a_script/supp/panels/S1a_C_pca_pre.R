#!/usr/bin/env Rscript
# S1a Figure C: PCA before normalization.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

pca_pre_df <- int_norm$pca_pre$scores |>
  left_join(int_norm$outlier_diag |> select(Col_ID, prefix),
    by = join_by(Col_ID)
  ) |>
  mutate(age = factor(ifelse(grepl("^O", prefix), "Old", "Young"),
    levels = c("Young", "Old")
  ))

pC <- ggplot(pca_pre_df, aes(PC1, PC2, color = age, fill = age)) +
  stat_ellipse(aes(group = age),
    geom = "polygon",
    alpha = 0.12, level = 0.80, linewidth = 0.3, show.legend = FALSE
  ) +
  geom_point(aes(shape = Timepoint), size = 1.6, alpha = 0.85) +
  scale_color_manual(values = AGE_COLORS, name = NULL) +
  scale_fill_manual(values = AGE_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, name = NULL) +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_pre$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_pre$var_exp[2]),
    tag = "C", title = "PCA before normalization",
    subtitle = sprintf(
      "%d samples | %.1f%% of variance in PC1-2 | 80%% ellipses",
      nrow(pca_pre_df), sum(int_norm$pca_pre$var_exp[1:2])
    )
  ) +
  FIG_THEME +
  theme(legend.position = "top", legend.key.size = unit(2.5, "mm"))

save_panel(pC, "S1a_C_pca_pre", "panel_C", pca_pre_df)

invisible(pC)
