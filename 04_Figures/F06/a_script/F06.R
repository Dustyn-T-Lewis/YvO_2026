#!/usr/bin/env Rscript
# Figure 6: A per-module ROCs for what age classifies, B per-module ROCs for
# what training classifies, younger over older, C module-phenotype coupling,
# six hero screen cells in a 3x2 scatter grid.

setwd(here::here())

pacman::p_load(patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F06/a_script/panels"
# A runs S7_A_module_grid.R first when the module grid is missing.
pA <- source_panel(file.path(PANELS, "A_age_roc.R")) +
  plot_annotation(title = NULL, subtitle = NULL)
pB <- source_panel(file.path(PANELS, "B_training_roc.R")) +
  plot_annotation(title = NULL, subtitle = NULL)
panel_C <- source_panel(file.path(PANELS, "C_hero_grid.R"))
pC_title <- panel_C$patches$annotation$title
pC_subtitle <- panel_C$patches$annotation$subtitle
pC <- panel_C +
  plot_annotation(title = NULL, subtitle = NULL)

# The key serves both panels, so it spans the canvas rather than sitting under
# one of them -- inside a half-width panel its labels were truncated.
pKey <- local({
  source(file.path(PANELS, "_roc_grid.R"), local = TRUE)
  key_composite
})

# The standalone A and B renders carry no title, so their titles are set here.
# Every cell in both panels already exists in module_grid_summary; the split
# changes which are drawn, not what was computed. Two cells that used to reach
# the top-twelve slice belong to neither question -- Pre vs Post pooled across
# ages, and the DeltaVL responder split -- so the subtitle says where they went
# rather than letting a BH-significant cell disappear without comment.
# One line each, measured against the 93 mm a half-canvas panel allows. The
# panel C subtitle clipped at 99% of its budget, so these sit near 80%.
pA_title <- "Modules that separate younger from older"
pA_subtitle <- "Mean Pre/Post ME | q<.05 solid, p<.05 dashed"

pB_title <- "Modules that separate pre from post training"
# Which row is which is not stated: every cell already says "Blue . Younger"
# or "Blue . Older", so a line repeating it is ink for nothing.
pB_subtitle <- "No older module reaches p<.05 | rest in S7 Figure"

RPT <- "04_Figures/F06/b_reports"

# Canvas, in millimetres throughout. The previous version derived every title
# position from a chain of fractions and fine adjustments; with a third panel
# to place, stating each edge in mm and converting once is easier to check.
COMP_W <- 210
AB_H <- 63 # the row holding panels A and B, side by side
KEY_H <- 14 # the module colour key, shared by A and B, two rows
C_H <- 88 # the scatter grid, full width
TAG_AB <- 11 # spacer above the A/B row: one subtitle line each
TAG_C <- 13 # spacer above panel C
COMP_H <- TAG_AB + AB_H + KEY_H + TAG_C + C_H

composite <- (
  plot_spacer() /
    ((pA | pB) + plot_layout(widths = c(1, 1))) /
    pKey /
    plot_spacer() /
    pC
) +
  plot_layout(heights = c(TAG_AB, AB_H, KEY_H, TAG_C, C_H))

txt <- composite_text_sizes(COMP_W, 178)

# The figure embeds at 165.1 x 147 mm, so the composite is read at about 79% of
# the size it is drawn. Titles are scaled up to survive that reduction.
TITLE_SCALE <- 1.42
TAG_SZ <- round(txt$tag * TITLE_SCALE, 1)
TTL_SZ <- round(txt$title * TITLE_SCALE, 1)
SUB_SZ <- round(txt$subtitle * TITLE_SCALE, 1)

# y is measured from the bottom, so a panel's title sits at the top of the
# spacer above it. Subtitles drop one title-line below.
mm_y <- function(mm) mm / COMP_H
mm_x <- function(mm) mm / COMP_W
SUB_DROP <- mm_y((TTL_SZ / 72) * 25.4 * 1.35)

TAG_X <- mm_x(4)
X_TTL <- mm_x(8.4)
TAG_Y_AB <- 1 - mm_y(1.5)
TAG_Y_C <- mm_y(C_H + TAG_C) - mm_y(1.5)
# Panel B's title starts at the midline, where its half of the row begins.
TAG_X_B <- 0.5 + mm_x(1)

composite <- composite & theme(plot.margin = margin(0, 4, 0, 5))

panel_head <- function(g, tag, title, subtitle, x, y) {
  g +
    draw_label(tag,
      x = x, y = y + mm_y(0.4), size = TAG_SZ, fontface = "bold",
      hjust = 0, vjust = 1
    ) +
    draw_label(title,
      x = x + X_TTL, y = y,
      size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1
    ) +
    # draw_label defaults to lineheight 0.9, which leaves two subtitle lines
    # overlapping in the descender band. 1.0 clears them.
    draw_label(subtitle,
      x = x + X_TTL, y = y - SUB_DROP,
      size = SUB_SZ, fontface = "bold.italic", colour = "grey40",
      hjust = 0, vjust = 1, lineheight = 1.0
    )
}

composite_final <- ggdraw(composite) |>
  panel_head("A", pA_title, pA_subtitle, TAG_X, TAG_Y_AB) |>
  panel_head("B", pB_title, pB_subtitle, TAG_X_B, TAG_Y_AB) |>
  panel_head("C", pC_title, pC_subtitle, TAG_X, TAG_Y_C)

graphics.off()
pdf_device <- get_pdf_device()

ggsave(file.path(RPT, "F06.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)
ggsave(file.path(RPT, "F06.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm",
  dpi = 300, limitsize = FALSE
)
message("F06 composite saved: [A: age ROC | B: training ROC] / [C: hero grid]")
