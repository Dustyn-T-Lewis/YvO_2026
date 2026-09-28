#!/usr/bin/env Rscript
# Figure 1: training volume, DXA lean body mass and VL thickness, laid out at
# double-column width (F01) and at single-column width (F01_single_col).

setwd(here::here())

pacman::p_load(ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F01/a_script/main/panels"
pA <- source_panel(file.path(PANELS, "A_training_volume.R"))
pB <- source_panel(file.path(PANELS, "B_dxa_lbm.R"))
pC <- source_panel(file.path(PANELS, "C_vl_thickness.R"))

# Titles and subtitles come off the panels and are redrawn on the canvas.
heads <- list(
  A = list(
    title = pA$labels$title,
    sub = paste0('italic("', pA$labels$subtitle, '")')
  ),
  B = list(title = pB[[1]]$labels$title, sub = attr(pB, "subtitle_expr")),
  C = list(title = pC[[1]]$labels$title, sub = attr(pC, "subtitle_expr"))
)
pA <- strip_for_composite(pA)
pair <- function(p) {
  (strip_for_composite(p[[1]]) | strip_for_composite(p[[2]])) +
    plot_layout(widths = c(0.65, 0.35))
}
pB_comp <- pair(pB)
pC_comp <- pair(pC)

# Draws one panel's tag, title and subtitle; annotate() for the subtitle
# because it is plotmath and draw_label() cannot parse an expression.
add_head <- function(canvas, tag, x, ttl_x, y, sub_off, sz) {
  canvas +
    draw_label(tag,
      x = x, y = y, size = sz$tag, fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(heads[[tag]]$title,
      x = ttl_x, y = y, size = sz$title, fontface = "bold", hjust = 0, vjust = 1
    ) +
    annotate("text",
      x = ttl_x, y = y - sub_off, label = heads[[tag]]$sub,
      parse = TRUE, hjust = 0, vjust = 1, size = sz$subtitle / .pt,
      colour = "grey30"
    )
}

RPT <- "04_Figures/F01/b_reports/main"

# Single-column layout, 85 x 125 mm; y holds each panel's top, normalized.
sc_cfg <- list(
  w = 85, h = 125, tag_x = 0.02, ttl_x = 0.08, sub_off = 0.022,
  y = c(A = 0.979, B = 0.614, C = 0.314)
)

# Double-column layout, 178 x 75 mm. Tags sit at x_a (A) and x_bc (B, C).
dc_cfg <- list(
  w = 178, h = 75,
  sz = list(
    tag = BASE_TAG, title = FIG_TITLE_SIZE, subtitle = FIG_SUBTITLE_SIZE
  ),
  x_a = 0.010, x_bc = 0.360,
  ttl_off = 0.030, # title x offset from tag
  sub_off = 0.026,
  # A and B sit on this row; C hangs off y_mid, so raising it lifts the two
  # panels asked for without moving C. The title-to-subtitle gap itself is
  # at its floor: at sub_off 0.024 the glyphs clear by 0.06 pt.
  y_top = 0.989,
  y_mid = 0.516
)

sc <- (pA / pB_comp / pC_comp) +
  plot_layout(heights = c(1.0, 0.8, 0.8)) &
  theme(plot.margin = margin(10, 2, 2, 2))

txt <- composite_text_sizes(sc_cfg$w)

sc <- ggdraw(sc)
for (k in c("A", "B", "C")) {
  sc <- add_head(
    sc, k, sc_cfg$tag_x, sc_cfg$ttl_x, sc_cfg$y[k], sc_cfg$sub_off, txt
  )
}

ggsave(file.path(RPT, "F01_single_col.pdf"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT, "F01_single_col.png"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", dpi = 300
)

dc <- wrap_elements(full = pA + theme(plot.margin = margin(8, 2, 10, 2))) +
  wrap_elements(full = pB_comp & theme(plot.margin = margin(4, 2, 2, 2))) +
  wrap_elements(full = pC_comp & theme(plot.margin = margin(4, 2, 4, 2))) +
  plot_layout(design = "AB\nAC", widths = c(0.35, 0.65), heights = c(1, 1))

dc_head <- function(canvas, tag, x, y) {
  add_head(canvas, tag, x, x + dc_cfg$ttl_off, y, dc_cfg$sub_off, dc_cfg$sz)
}
dc <- ggdraw(dc) |>
  dc_head("A", dc_cfg$x_a, dc_cfg$y_top) |>
  dc_head("B", dc_cfg$x_bc, dc_cfg$y_top) |>
  dc_head("C", dc_cfg$x_bc, dc_cfg$y_mid)

ggsave(file.path(RPT, "F01.pdf"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT, "F01.png"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", dpi = 300
)

message("F01 main composites done")
