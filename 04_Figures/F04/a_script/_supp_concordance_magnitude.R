#!/usr/bin/env Rscript
# SUPP Panel F: sign concordance against effect magnitude
# Diagnostic for main Panel A — asks whether whole-proteome sign discordance
# carries biology or just tracks proteins too small to have a reliable sign.
# Proteins are binned on the smaller of the two age-group |log2FC| values, so a
# bin is low only when the effect is near zero in both groups.
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Exports: pS_conc_mag (ggplot)

BASE <- "04_Figures/F04"
DAT <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

conc_df <- dep |>
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) |>
  mutate(
    concordant = sign(logFC_Training_Young) == sign(logFC_Training_Old),
    min_abs = pmin(abs(logFC_Training_Young), abs(logFC_Training_Old)),
    young_fdr_dep = !is.na(adj.P.Val_Training_Young) &
      adj.P.Val_Training_Young < 0.05
  )

quint <- quantile(conc_df$min_abs, seq(0, 1, 0.2))
conc_df <- conc_df |>
  mutate(bin = cut(min_abs, quint, include.lowest = TRUE, labels = FALSE))

bin_df <- conc_df |>
  group_by(bin) |>
  summarise(
    n = n(),
    lo = min(min_abs), hi = max(min_abs),
    n_conc = sum(concordant),
    pct = 100 * mean(concordant),
    .groups = "drop"
  ) |>
  mutate(
    label = sprintf("%.2f-%.2f", lo, hi),
    se = 100 * sqrt((pct / 100) * (1 - pct / 100) / n)
  )

dep_young_fdr <- filter(conc_df, young_fdr_dep)
pct_young_fdr <- 100 * mean(dep_young_fdr$concordant)
trend_p <- prop.trend.test(bin_df$n_conc, bin_df$n)$p.value

pS_conc_mag <- ggplot(bin_df, aes(x = factor(bin), y = pct)) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey55") +
  geom_hline(
    yintercept = pct_young_fdr, linetype = "dotted",
    colour = AGE_COLORS[["Old"]], linewidth = 0.5
  ) +
  geom_col(fill = AGE_COLORS[["Young"]], width = 0.7) +
  geom_errorbar(aes(ymin = pct - se, ymax = pct + se),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_text(aes(y = pct + se, label = sprintf("%.0f%%", pct)),
    vjust = -0.8, size = 1.7
  ) +
  annotate("text",
    x = 0.6, y = 51.5, hjust = 0, size = 1.7, colour = "grey40",
    label = "chance"
  ) +
  annotate("text",
    x = 0.6, y = pct_young_fdr + 1.5, hjust = 0, size = 1.7,
    colour = AGE_COLORS[["Old"]],
    label = sprintf("%d FDR DEPs (%.0f%%)", nrow(dep_young_fdr), pct_young_fdr)
  ) +
  scale_x_discrete(labels = bin_df$label) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Sign concordance by effect magnitude",
    subtitle = sprintf(
      "Quintiles of min(|log2FC|) across age groups; trend %s",
      fmt_p(trend_p)
    ),
    x = "min(|log2FC|) across age groups", y = "Concordant (%)"
  ) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))

write_csv(bin_df, file.path(DAT, "SUPP_concordance_magnitude.csv"))

RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89
PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_concordance_magnitude.png"), pS_conc_mag,
  width = PW, height = PH, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_concordance_magnitude.pdf"), pS_conc_mag,
  width = PW, height = PH, units = "mm", device = pdf_device
)

pS_conc_mag <- strip_for_composite(pS_conc_mag)

message(sprintf(
  "SUPP concordance-magnitude: bottom quintile %.1f%%, top quintile %.1f%%, %d FDR DEPs %.1f%%, trend p = %.3g",
  bin_df$pct[1], bin_df$pct[nrow(bin_df)], nrow(dep_young_fdr), pct_young_fdr, trend_p
))
