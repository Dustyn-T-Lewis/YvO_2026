#!/usr/bin/env Rscript
# Figure 2: PCA, log2FC distributions, DEP counts, UpSet overlap, fGSEA counts
# and the t-statistic barcode in a 3 x 2 grid.

setwd(here::here())

pacman::p_load(ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F02/a_script/panels"
pA <- source_panel(file.path(PANELS, "A_pca.R"))
pB <- source_panel(file.path(PANELS, "B_logfc_density.R"))
pC <- source_panel(file.path(PANELS, "C_dep_counts.R"))
pD_standalone <- source_panel(file.path(PANELS, "D_upset.R"))
pE <- source_panel(file.path(PANELS, "E_fgsea.R"))
pF <- source_panel(file.path(PANELS, "F_barcode.R"))

pA_title <- pA$labels$title
pA_subtitle <- pA$labels$subtitle
pA <- pA + labs(title = NULL, subtitle = NULL, tag = NULL)

pB_title <- pB$labels$title
pB_subtitle <- pB$labels$subtitle %||% ""
pB <- strip_for_composite(pB)

pC_title <- pC[[1]]$labels$title
pC_subtitle <- pC[[1]]$labels$subtitle
pC <- pC & labs(title = NULL, subtitle = NULL, tag = NULL) & theme(legend.position = "none")

pD_bars <- attr(pD_standalone, "bars")
pD_dots <- attr(pD_standalone, "dots")
p_key_dir_D <- attr(pD_standalone, "key")
pD_title <- pD_bars$labels$title
pD_subtitle <- pD_bars$labels$subtitle
# Build clean version (no title/subtitle) for composite export.
#
# The composite copy also drops the bar plot's 5 pt left margin and trims the
# patchwork's own left margin: together they pull the shared left edge 8.25 pt
# out so D's panel border lines up with panel A's at ~22.3 pt. Only the
# composite needs this -- standalone D keeps the margin that stops "100" from
# clipping against the device edge. "Intersection size" is drawn on the
# canvas, not the gtable, so it stays where it is; the bar plot's own y-axis
# numbers move left to within ~2 pt of it, which is as far as this can go.
pD_bars_clean <- pD_bars + labs(title = NULL, subtitle = NULL) +
  theme(plot.margin = margin(2, 0, 0, 0))
pD_pw <- (pD_bars_clean / pD_dots) + plot_layout(heights = c(0.78, 0.22)) +
  plot_annotation(
    theme = theme(plot.margin = margin(t = 4, r = 2, b = 4, l = 1.75))
  )

# Clean version without title or subtitle; the composite places the title
pD <- ggdraw(pD_pw) +
  draw_label("Intersection size",
    x = 0.02, y = 0.58, angle = 90,
    size = 5, fontface = "bold"
  ) +
  # Drawn on the whole-panel canvas, so y has to clear the dot matrix and the
  # bar tops. 0.757 puts its first label level with panel E's two key columns:
  # E's sit at 171.8 pt on the composite, D's at 181.6, and the bottom row is
  # 172.77 pt tall, so 9.8 / 172.77 = 0.057 above the old 0.70.
  draw_plot(p_key_dir_D, x = 0.87, y = 0.757, width = 0.12, height = 0.22)

pE_title <- pE[[1]]$labels$title
pE_subtitle <- pE[[1]]$labels$subtitle
# pE is a patchwork with inset_element legends, so & strips all its plots.
# Inset legends use theme_void() so stripping is harmless to them.
pE <- pE & labs(title = NULL, subtitle = NULL, tag = NULL)

pF_title <- pF$labels$title
pF_subtitle <- pF$labels$subtitle
pF <- strip_for_composite(pF)

layout <- "ABC\n###\nDEF"
# Top row deliberately deeper than half so the D/E/F titles and subtitles get
# clear air below the top row's axis labels. BOT_Y below moves with it.
ROW_TOP <- 0.470
SPACER <- 0.00

# Per-panel margins (top breathing room + per-panel width nudges).
# Bottom margin matched to pB/pC's 0-8pt range -- the previous 12pt left A's
# rendered plot visibly shorter than B's in the composite row.
pA <- pA + theme(plot.margin = margin(12, 2, 7, 2))
pB <- pB + theme(plot.margin = margin(12, -52, 0, 5))
# Index into just the base bar plot -- `&` broadcasts to every patch,
# including the key inset, and overwrote its own zero-margin theme, pushing
# the key down onto the Tr.(O)/Interaction row boundary.
pC[[1]] <- pC[[1]] + theme(plot.margin = margin(6, 2, 3, 0))
# Top margin no longer negative: the deeper ROW_TOP shortens the bottom row,
# and the old -3 pulled F's first callout up into its own subtitle.
pF <- pF + theme(plot.margin = margin(9, 2, -12, 0))
# D and E started their plots 12 and 14 px above F's top border, close enough
# to their subtitles to read as crowded. Only the top moves, so both lose that
# much height. pE indexes like pC: `&` would reach its inset legends too.
pD <- pD + theme(plot.margin = margin(3, 0, 0, 0))
pE[[1]] <- pE[[1]] + theme(plot.margin = margin(3.5, 0, 0, 0))

composite <- wrap_elements(full = pA) + pB + wrap_elements(full = pC) +
  wrap_elements(full = pD) +
  wrap_elements(full = pE) +
  pF +
  plot_layout(
    design  = layout,
    widths  = c(160, 127, 138),
    heights = c(ROW_TOP, SPACER, 1 - ROW_TOP - SPACER)
  )

# Manual tag + title + subtitle placement via cowplot.
COMP_W <- 178
COMP_H <- 115
TAG_SZ <- composite_text_sizes(COMP_W)$tag
TTL_SZ <- composite_text_sizes(COMP_W)$title
SUB_SZ <- composite_text_sizes(COMP_W)$subtitle
TOP_Y <- 0.995 - 2 / COMP_H + 0.020 + 0.002 - 0.005 - 0.009
# Baseline of the D/E/F title row; tracks ROW_TOP, so shift both by the same
# amount whenever the bottom row moves.
BOT_Y <- 0.528
X_LEFT <- 0.002
X_MID <- 0.372
X_RIGHT <- 0.630
X_MID_BOT <- X_MID
X_RIGHT_BOT <- X_RIGHT + 0.010
X_TTL <- 0.04
TTL_NUDGE <- -0.008
BE_NUDGE <- 0.021
TAG_DY <- -0.002
SUB_OFFSET <- 0.022

composite <- ggdraw(composite) +
  # Panel A
  draw_label("A", x = X_LEFT, y = TOP_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = X_LEFT + X_TTL, y = TOP_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = X_LEFT + X_TTL, y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel B
  draw_label("B", x = X_MID, y = TOP_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = X_MID + X_TTL - BE_NUDGE, y = TOP_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = X_MID + X_TTL - BE_NUDGE, y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel C
  draw_label("C", x = X_RIGHT + 0.042, y = TOP_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = TOP_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel D
  draw_label("D", x = X_LEFT, y = BOT_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pD_title, x = X_LEFT + X_TTL, y = BOT_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pD_subtitle, x = X_LEFT + X_TTL, y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel E
  draw_label("E", x = X_MID_BOT, y = BOT_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pE_title, x = X_MID + X_TTL - BE_NUDGE, y = BOT_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pE_subtitle, x = X_MID + X_TTL - BE_NUDGE, y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel F
  draw_label("F", x = X_RIGHT_BOT + 0.042, y = BOT_Y - TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pF_title, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = BOT_Y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pF_subtitle, x = X_RIGHT + 0.042 + X_TTL + TTL_NUDGE, y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

RPT <- "04_Figures/F02/b_reports"
pdf_dev <- get_pdf_device()

ggsave(file.path(RPT, "F02.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev
)
ggsave(file.path(RPT, "F02.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message("F02 main composite done")
