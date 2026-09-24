#!/usr/bin/env Rscript
# Figure 2A: sample PCA with bootstrap CIs on the variance explained and the
# PERMANOVA label.

setwd(here::here())
source("04_Figures/F02/a_script/main/panels/_main.R", local = TRUE)
PC_W <- 67
PC_H <- 55

pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
boot_var <- replicate(1000, {
  idx <- sample(nrow(imp_mat), replace = TRUE)
  100 * summary(prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE))$importance[2, 1:2]
})
var_ci <- data.frame(
  PC = c("PC1", "PC2"), var_pct = var_pct,
  ci_lo = apply(boot_var, 1, quantile, 0.025),
  ci_hi = apply(boot_var, 1, quantile, 0.975)
)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(imp_meta, by = join_by(sample_id))

dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ age * time,
  data = imp_meta,
  permutations = how(nperm = 999, blocks = imp_meta$subject),
  by = "terms"
)
perm_terms <- c("age", "time", "age:time")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  " PERMANOVA\nAge        R\u00b2 = %.3f,  %s\nTime       R\u00b2 = %.3f,  %s\nAge\u00d7Time R\u00b2 = %.3f,  %s",
  perm_r2[1], fmt_p(perm_pv[1]), perm_r2[2], fmt_p(perm_pv[2]),
  perm_r2[3], fmt_p(perm_pv[3])
)

bd_age_p <- permutest(betadisper(dist_mat, imp_meta$age), permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(betadisper(dist_mat, imp_meta$time), permutations = 999)$tab$`Pr(>F)`[1]
if (bd_age_p < 0.05 || bd_time_p < 0.05) {
  warning("Heterogeneous dispersions — interpret PERMANOVA with caution")
}

pA <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group),
    geom = "polygon",
    alpha = 0.10, level = 0.80, show.legend = FALSE
  ) +
  stat_ellipse(aes(group = group),
    level = 0.80, linewidth = 0.4,
    linetype = "dashed", show.legend = FALSE
  ) +
  geom_point(size = 1.8, alpha = 0.85) +
  annotate("label",
    x = -Inf, y = Inf, label = perm_label,
    hjust = -0.02, vjust = 1.05, lineheight = 0.9,
    size = STAT_BOX_PT / .pt, color = "grey20", fontface = "bold.italic",
    fill = alpha("white", 0.85), linewidth = 0.2,
    label.padding = unit(0.10, "lines")
  ) +
  scale_color_manual(
    values = PCA_COLORS,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post"),
    guide = guide_legend(override.aes = list(size = 1.6))
  ) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(
    values = PCA_SHAPES,
    labels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
  ) +
  labs(
    title = "Sample PCA",
    subtitle = sprintf(
      "n = %d, %s proteins (imputed)",
      nrow(imp_meta), format(nrow(imp_mat), big.mark = ",")
    ),
    x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
    y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2])
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE,
      face = "bold.italic", color = "grey30"
    ),
    legend.position = c(0.86, 0.14), legend.background = element_blank(),
    legend.key = element_blank(), legend.title = element_blank(),
    # Matched to panel C's hand-rolled key (geom_text size 1.25 mm = 3.56 pt)
    # so the two group keys carry the same visual weight across the row.
    legend.text = element_text(size = 3.6, face = "bold"),
    legend.key.size = unit(1.8, "mm"),
    plot.margin = margin(6, 6, 2, 8)
  )

write.csv(var_ci, file.path(DAT, "panel_A_pca_variance_ci.csv"), row.names = FALSE)
ggsave(file.path(PNL_PNG, "A_pca.png"), pA,
  width = PC_W, height = PC_H, units = "mm", dpi = 300
)
ggsave(file.path(PNL_PDF, "A_pca.pdf"), pA,
  width = PC_W, height = PC_H, units = "mm", device = pdf_dev
)

invisible(pA)
