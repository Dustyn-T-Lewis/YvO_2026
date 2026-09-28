#!/usr/bin/env Rscript
# Figure 3: the four volcano rings in a 2 x 2 grid with the NES and opacity keys.

setwd(here::here())

pacman::p_load(readr, dplyr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

PANELS <- "04_Figures/F03/a_script/main/panels"
pA <- strip_for_composite(source_panel(file.path(PANELS, "A_volcano_aging.R")))
pB <- strip_for_composite(source_panel(file.path(PANELS, "B_volcano_young.R")))
pC <- strip_for_composite(source_panel(file.path(PANELS, "C_volcano_old.R")))
pD <- strip_for_composite(source_panel(file.path(PANELS, "D_volcano_interaction.R")))

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
fgsea_all <- read_csv("04_Figures/shared/fgsea_tstat_all_v2.csv", show_col_types = FALSE)

PI_THRESH <- 0.05
FDR_THRESH <- 0.05

RPT <- "04_Figures/F03/b_reports/main"

# The ring draws on Hallmark + GO Slim (select_ring_terms() default), so the
# subtitle counts the same two databases rather than the full five.
RING_DBS <- c("Hallmark", "GO Slim")
# Pi counts are absolute, not Pi-only: a protein clearing both criteria is in
# both totals. The faded points are the Pi-only subset, which the key at the
# foot of the figure names rather than the subtitle.
contrast_stats <- function(ctr) {
  fdr <- dep_df[[paste0("adj.P.Val_", ctr)]]
  pi <- dep_df[[paste0("pi_score_", ctr)]]
  lfc <- dep_df[[paste0("logFC_", ctr)]]
  ok <- !is.na(fdr) & !is.na(pi) & !is.na(lfc)
  ctr_rows <- fgsea_all[fgsea_all$contrast == ctr & fgsea_all$database %in% RING_DBS, ]
  # Two lines: at 9.5 pt a half-width column holds about 44 characters, and the
  # one-line form ran into the neighbouring panel's subtitle.
  sprintf(
    "FDR %d up, %d down  |  \u03a0 %d up, %d down\n%d / %d pathways enriched",
    sum(ok & fdr < FDR_THRESH & lfc > 0), sum(ok & fdr < FDR_THRESH & lfc <= 0),
    sum(ok & pi < PI_THRESH & lfc > 0), sum(ok & pi < PI_THRESH & lfc <= 0),
    sum(ctr_rows$padj < FDR_THRESH, na.rm = TRUE), sum(!is.na(ctr_rows$padj))
  )
}

# The two point opacities the rings use, named. point_alpha is 0.55 in
# _build_volcano_panel.R and volcano_ring.R scales it by 1.4 and 0.45.
tier_key <- local({
  tiers <- c(min(0.55 * 1.4, 1), 0.55 * 0.45)
  kdf <- data.frame(
    x = rep(c(0, 0.30), each = 2),
    y = rep(c(1, 0), 2),
    a = rep(tiers, 2),
    col = rep(c(DIR_COLORS[["Up"]], DIR_COLORS[["Down"]]), each = 2)
  )
  ldf <- data.frame(y = c(1, 0), lab = c("FDR < 0.05", "\u03a0 < 0.05 only"))
  ggplot(kdf, aes(x = x, y = y)) +
    geom_point(aes(alpha = a, colour = col), size = 2.4) +
    geom_text(data = ldf, aes(y = y, label = lab), x = 0.52, hjust = 0,
              size = 7 / .pt, fontface = "bold", colour = "grey25",
              inherit.aes = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_x_continuous(limits = c(-0.2, 3.1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.7, 1.7), expand = c(0, 0)) +
    theme_void()
})

nes_legend <- build_nes_legend_bar(text_size = 6.5, title_size = 6.5,
                                   bar_margin = margin(0, 0, 0, 0, "mm"))

# The title block is ~10 pt taller per row than the old 7/6 pt one, so every
# ring drops by that much. Panel B binds: its "Protein Folding" label sits
# closest to 12 o'clock and grazed the old subtitle.
pA <- pA + theme(plot.margin = margin(13, -9, 0, 9, "mm"))
pB <- pB + theme(plot.margin = margin(13, 0, 0, 0, "mm"))
pC <- pC + theme(plot.margin = margin(8, -9, 0, 9, "mm"))
pD <- pD + theme(plot.margin = margin(8, 0, 0, 0, "mm"))

composite <- ((pA | pB) / (pC | pD)) + plot_layout(heights = c(1, 1))

# 3 mm of extra height pays back most of the ring area the deeper title block
# takes, without pushing the embedded figure past its 165.1 mm column width.
COMP_W <- 178; COMP_H <- 190
# Set explicitly rather than via composite_text_sizes(COMP_W). That helper keys
# off canvas *width*, which is right for the wide, short canvases of the other
# figures but not here: F03 is near-square, so at the shared 165.1 mm embed width
# it stands 169.7 mm tall and its 7 pt title is only 7 / 518.7 = 1.4% of figure
# height, against 2.0-3.3% in Figures 1, 2, 4 and 5. Hence the apparent shrink.
# 11 pt over the 183 mm (518.7 pt) canvas is 2.12%, back inside that band, and
# the subtitle and tag keep their 6:7 and 8:7 ratios to it.
TAG_SZ <- 12
TTL_SZ <- 11
SUB_SZ <- 9.5
# Tags for the left (A, C) and right (B, D) columns; titles sit X_TTL right of
# the tag and subtitles SUB_OFF below the title, tight but clear of the rings.
X_L <- 0.070; X_R <- 0.510; X_TTL <- 0.040
Y_TOP <- 0.9647; Y_BOT <- 0.4947; SUB_OFF <- 0.0247

titles <- c("Aging Effect", "Training Response (Young)",
            "Training Response (Old)", "Age \u00d7 Training Interaction")
subs <- vapply(c("Aging", "Training_Young", "Training_Old", "Interaction"),
               contrast_stats, character(1), USE.NAMES = FALSE)
tags <- LETTERS[1:4]
xs <- c(X_L, X_R, X_L, X_R)
ys <- c(Y_TOP, Y_TOP, Y_BOT, Y_BOT)

composite <- ggdraw(composite)
for (i in 1:4) {
  composite <- composite +
    draw_label(tags[i], x = xs[i], y = ys[i] + 0.002, size = TAG_SZ,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(titles[i], x = xs[i] + X_TTL, y = ys[i], size = TTL_SZ,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(subs[i], x = xs[i] + X_TTL, y = ys[i] - SUB_OFF, size = SUB_SZ,
               fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1)
}
composite <- composite +
  draw_plot(nes_legend, x = 0.22, y = 0.022, width = 0.30, height = 0.034) +
  draw_plot(tier_key, x = 0.55, y = 0.018, width = 0.26, height = 0.044)

ggsave(file.path(RPT, "F03.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device(),
       limitsize = FALSE)
ggsave(file.path(RPT, "F03.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)

message("F03 main composite done")
