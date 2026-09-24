#!/usr/bin/env Rscript
# Figure 2B: log2 fold-change distributions with the KS and blunting statistics.

setwd(here::here())
source("04_Figures/F02/a_script/main/panels/_main.R", local = TRUE)
PD_W <- 48
PD_H <- 55

lfc_long_all <- dep_df |>
  select(any_of(c("uniprot_id", "gene")), starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC))
write_csv(lfc_long_all, file.path(DAT, "panel_B_logfc_long.csv"))

lfc_long <- lfc_long_all |>
  filter(contrast %in% c("Aging", "Training_Young", "Training_Old")) |>
  mutate(contrast = factor(contrast, levels = c("Aging", "Training_Young", "Training_Old")))

set.seed(42)
lfc_stats <- lfc_long |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    boot_median_ci(abs(logFC)),
    n_above_05 = sum(abs(logFC) > 0.5), .by = contrast
  )
write_csv(lfc_stats, file.path(DAT, "panel_B_stats.csv"))

lfc_stats$annotation <- sprintf(
  "Med.|logFC| = %.2f\nn(>0.5) = %d",
  lfc_stats$med_abs_lfc, lfc_stats$n_above_05
)

dist_subtitle <- NULL
blunt <- tryCatch(as.data.frame(read_excel(DEP_XLSX, sheet = "blunting")), error = \(e) NULL)
if (!is.null(blunt)) {
  ks_row <- blunt[grepl("KS", blunt$test), ]
  # Resampled by protein, not by column. The two contrasts are paired on the
  # same rows and correlated at rho ~ 0.3; drawing them independently breaks
  # that pairing and misstates the interval. Pivoting first also guarantees the
  # two vectors are gene-aligned, which filtering the long frame does not.
  blunt_pair <- lfc_long |>
    filter(contrast %in% c("Training_Young", "Training_Old")) |>
    pivot_wider(
      id_cols = any_of(c("uniprot_id", "gene")),
      names_from = contrast, values_from = logFC
    ) |>
    filter(!is.na(Training_Young), !is.na(Training_Old)) |>
    transmute(y = abs(Training_Young), o = abs(Training_Old))
  obs_ratio <- median(blunt_pair$o) / median(blunt_pair$y)
  set.seed(42)
  ratio_ci <- quantile(
    replicate(2000, {
      i <- sample(nrow(blunt_pair), replace = TRUE)
      median(blunt_pair$o[i]) / median(blunt_pair$y[i])
    }),
    c(0.025, 0.975)
  )
  dist_subtitle <- sprintf(
    "KS D = %.2f, %s | blunting %.2f [%.2f, %.2f]",
    ks_row$statistic, fmt_p(ks_row$p_value),
    obs_ratio, ratio_ci[1], ratio_ci[2]
  )
}

lfc_binwidth <- 4 / 50
pB <- ggplot(lfc_long, aes(logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
    alpha = 0.15, linewidth = 0.5, color = "grey20"
  ) +
  geom_text(
    data = lfc_stats,
    aes(x = 0, y = 695, label = CTR_SHORT[as.character(contrast)]),
    inherit.aes = FALSE, hjust = 0.5, vjust = 1,
    size = FIG_LEGEND_TEXT / .pt, color = "grey30", fontface = "bold"
  ) +
  geom_label(
    data = lfc_stats,
    aes(x = -1.08, y = 695, label = annotation),
    inherit.aes = FALSE, hjust = 0, vjust = 1,
    size = STAT_BOX_PT / .pt, color = "grey20", fontface = "bold.italic",
    lineheight = 0.9, fill = alpha("white", 0.85), linewidth = 0.2,
    label.padding = unit(0.10, "lines")
  ) +
  facet_wrap(~contrast, ncol = 1, labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 700)) +
  scale_fill_manual(values = CONTRAST_COLORS[c("Aging", "Training_Young", "Training_Old")]) +
  labs(
    title = "Effect Size Distribution", subtitle = dist_subtitle,
    x = expression(bold(log[2] ~ FC)), y = " "
  ) +
  FIG_THEME +
  theme(
    legend.position = "none", strip.text = element_blank(),
    strip.background = element_blank(), panel.spacing.y = unit(0, "pt"),
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE,
      face = "bold.italic", color = "grey40"
    ),
    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
    axis.ticks.y = element_blank(),
    plot.margin = margin(6, 4, 0, 4)
  )
ggsave(file.path(PNL_PNG, "B_logfc_density.png"), pB,
  width = PD_W, height = PD_H, units = "mm", dpi = 300
)
ggsave(file.path(PNL_PDF, "B_logfc_density.pdf"), pB,
  width = PD_W, height = PD_H, units = "mm", device = pdf_dev
)

invisible(pB)
