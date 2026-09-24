#!/usr/bin/env Rscript
# S1a Figure D: PCA after cyclic loess.

setwd(here::here())
source("04_Figures/F00/a_script/supp/panels/_inputs.R", local = TRUE)

pD <- ggplot(
  int_norm$pca_post$scores,
  aes(PC1, PC2, color = Group_Time, fill = Group_Time)
) +
  stat_ellipse(aes(group = Group_Time),
    geom = "polygon",
    alpha = 0.12, level = 0.80, linewidth = 0.3, show.legend = FALSE
  ) +
  geom_point(aes(shape = Timepoint), size = 1.6, alpha = 0.85) +
  scale_color_manual(
    values = PCA_COLORS, name = NULL,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  labs(
    x = sprintf("PC1 (%.1f%%)", int_norm$pca_post$var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", int_norm$pca_post$var_exp[2]),
    tag = "D", title = "PCA after cyclic loess",
    subtitle = sprintf(
      "%d samples | %.1f%% of variance in PC1-2 | 80%% ellipses",
      nrow(int_norm$pca_post$scores), sum(int_norm$pca_post$var_exp[1:2])
    )
  ) +
  FIG_THEME +
  theme(
    legend.position = "top", legend.key.size = unit(2.5, "mm"),
    legend.text = element_text(size = 5.5)
  )

save_panel(pD, "S1a_D_pca_post", "panel_D", int_norm$pca_post$scores)

invisible(pD)
