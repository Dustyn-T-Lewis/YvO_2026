#!/usr/bin/env Rscript
# Figure 2C: DEPs per contrast as a percentage of the proteome, at FDR and Pi.

setwd(here::here())
source("04_Figures/F02/a_script/panels/_main.R", local = TRUE)
PA_W <- 67
PA_H <- 55
all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

sig_sets <- list()
dir_map <- list()
for (ctr in CONTRASTS) {
  pi_vals <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]] <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"), dep_df$gene[is_sig])
}
n_total <- length(all_genes)

pi_total <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, \(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.05, na.rm = TRUE) else 0
}))
SET_DISPLAY_COLORS <- c(
  "Aging" = unname(CONTRAST_COLORS["Aging"]),
  "Tr.(Y)" = unname(CONTRAST_COLORS["Training_Young"]),
  "Tr.(O)" = unname(CONTRAST_COLORS["Training_Old"]),
  "Inter." = unname(CONTRAST_COLORS["Interaction"])
)

frac_df <- bind_rows(lapply(CONTRASTS, \(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  tibble(
    contrast = SET_LABELS[ctr], threshold = c("q < 0.05", "Π < 0.05"),
    n = c(
      sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05),
      length(sig_sets[[ctr]])
    )
  )
})) |>
  mutate(
    contrast = factor(contrast, levels = rev(c("Aging", "Tr.(Y)", "Tr.(O)", "Inter."))),
    threshold = factor(threshold, levels = c("q < 0.05", "Π < 0.05")),
    pct = 100 * n / n_total, fill_key = paste(contrast, threshold, sep = "___")
  ) |>
  arrange(contrast, desc(pct)) |>
  filter(n > 1)

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "q < 0.05", sep = "___")] <- adjustcolor(col, alpha.f = 0.40)
  FRAC_FILL[paste(cname, "Π < 0.05", sep = "___")] <- col
}

# Axis fit to this engine's actual max rather than a fixed span -- a shared
# 0-28% range across engines left every non-Aging bar squeezed into a sliver.
# Headroom covers the stacked "n FDR / n \u03A0" label at the longer bar's end.
y_max_C <- max(frac_df$pct) * 1.22

# One stacked label per contrast anchored at whichever bar reaches further,
# rather than separately-positioned numbers. Built straight from
# dep_df/sig_sets, not frac_df, so a tier's count still shows even when its own
# bar is too small (n <= 1) to plot. The nominal-p tier was dropped for R1.3
# and R2.1a: reporting it alongside FDR invited the reader to treat it as a
# result, which is what both reviewers objected to.
label_df_C <- tibble(
  contrast = factor(SET_LABELS[CONTRASTS], levels = levels(frac_df$contrast)),
  n_fdr = sapply(CONTRASTS, \(ctr) {
    fdr_col <- paste0("adj.P.Val_", ctr)
    sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.05)
  }),
  n_pi = sapply(CONTRASTS, \(ctr) length(sig_sets[[ctr]]))
) |>
  mutate(
    x_pos = pmax(100 * n_fdr / n_total, 100 * n_pi / n_total),
    label = sprintf("%d FDR\n%d \u03A0", n_fdr, n_pi)
  )

# Shade key: light = FDR, solid = Pi. Generic grey
# squares -- the bars themselves already carry per-contrast hue via the
# background wash. Small and tucked into the bottom-right corner so it reads
# as a footnote, not a fourth thing competing with the bars for attention.
key_df_C <- tibble(
  label = c("FDR < 0.05", "\u03A0 < 0.05"),
  y = c(0, -0.13),
  fill = c(adjustcolor("grey20", alpha.f = 0.40), "grey20")
)
p_key_C <- ggplot(key_df_C) +
  geom_point(aes(x = 0, y = y),
    shape = 22, size = 1.5,
    fill = key_df_C$fill, color = "grey20", stroke = 0.3
  ) +
  geom_text(aes(x = 0.29, y = y, label = label),
    size = 1.25, color = "grey20", hjust = 0, fontface = "bold"
  ) +
  scale_x_continuous(limits = c(-0.15, 2.1)) +
  coord_cartesian(ylim = c(-0.33, 0.05), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

pC <- ggplot(frac_df, aes(contrast, pct, fill = fill_key)) +
  annotate("rect",
    xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
    fill = CONTRAST_COLORS["Aging"], alpha = 0.20, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
    fill = CONTRAST_COLORS["Training_Young"], alpha = 0.20, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
    fill = CONTRAST_COLORS["Training_Old"], alpha = 0.20, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
    fill = CONTRAST_COLORS["Interaction"], alpha = 0.20, color = "grey70", linewidth = 0.2
  ) +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  # Exact counts stacked as one three-line block ("x p" / "y FDR" / "z Π")
  # at whichever bar reaches further -- the key gives the shade meaning,
  # this gives the numbers, read together rather than as floating figures.
  geom_text(
    data = label_df_C,
    aes(x = contrast, y = x_pos + y_max_C * 0.012, label = label),
    hjust = 0, vjust = 0.5, size = 1.5, fontface = "bold", color = "grey15",
    lineheight = 0.85, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = scales::pretty_breaks(n = 5), limits = c(0, y_max_C)
  ) +
  # Tightened from the discrete default (add = 0.6) -- the extra padding sat
  # as blank space above/below the outermost bars, making the panel's
  # rendered content noticeably shorter than A/B's in the composite row.
  scale_x_discrete(expand = expansion(add = 0.3)) +
  coord_flip() +
  labs(
    title = "DEPs per Contrast",
    subtitle = sprintf(
      "%s genes | FDR %d | \u03A0 %d",
      format(n_total, big.mark = ","), fdr_total, pi_total
    ),
    x = NULL, y = "% of proteome"
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(
      size = FIG_SUBTITLE_SIZE,
      face = "bold.italic", color = "grey40"
    ),
    legend.position = "none",
    axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT - 0.5)
  )
pC <- (pC + inset_element(p_key_C,
  left = 0.80, right = 0.99, top = 0.20, bottom = 0.02
)) & theme(legend.position = "none")

ggsave(file.path(PNL_PNG, "C_dep_counts.png"), pC,
  width = PA_W, height = PA_H, units = "mm", dpi = 300
)
ggsave(file.path(PNL_PDF, "C_dep_counts.pdf"), pC,
  width = PA_W, height = PA_H, units = "mm", device = pdf_dev
)

invisible(pC)
