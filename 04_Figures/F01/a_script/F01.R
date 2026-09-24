#!/usr/bin/env Rscript
# Figure 1: training volume, DXA lean body mass and VL thickness, laid out at
# double-column width (F01) and at single-column width (F01_single_col).

setwd(here::here())

pacman::p_load(ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F01/a_script/panels"
pA <- source_panel(file.path(PANELS, "A_training_volume.R"))
pB <- source_panel(file.path(PANELS, "B_dxa_lbm.R"))
pC <- source_panel(file.path(PANELS, "C_vl_thickness.R"))

pA_title <- pA$labels$title
pA_subtitle <- paste0('italic("', pA$labels$subtitle, '")')
pA <- strip_for_composite(pA)
pB_title <- pB[[1]]$labels$title
pB_subtitle <- attr(pB, "subtitle_expr")
pB_left <- strip_for_composite(pB[[1]])
pB_right <- strip_for_composite(pB[[2]])
pC_title <- pC[[1]]$labels$title
pC_subtitle <- attr(pC, "subtitle_expr")
pC_left <- strip_for_composite(pC[[1]])
pC_right <- strip_for_composite(pC[[2]])

RPT <- "04_Figures/F01/b_reports"

# Single-column layout (85 × 125 mm) — for journal column width
sc_cfg <- list(
  w = 85, h = 125,
  tag_x = 0.02,
  ttl_x = 0.08,
  sub_off = 0.022,
  y = c(A = 0.979, B = 0.614, C = 0.314) # panel top positions (normalized)
)

# Double-column layout (178 × 75 mm) — for full-width display
dc_cfg <- list(
  w = 178, h = 75,
  tag_sz = BASE_TAG,
  ttl_sz = FIG_TITLE_SIZE,
  sub_sz = FIG_SUBTITLE_SIZE,
  x_a = 0.010, # tag x for panel A column
  x_bc = 0.360, # tag x for panels B/C column
  ttl_off = 0.030, # title x offset from tag
  sub_off = 0.026,
  # A and B sit on this row; C hangs off y_mid, so raising it lifts the two
  # panels asked for without moving C. The title-to-subtitle gap itself is
  # at its floor: at sub_off 0.024 the glyphs clear by 0.06 pt.
  y_top = 0.989,
  y_mid = 0.516 # mid row y (panel C)
)

pB_comp <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

sc <- (pA / pB_comp / pC_comp) +
  plot_layout(heights = c(1.0, 0.8, 0.8)) &
  theme(plot.margin = margin(10, 2, 2, 2))

txt <- composite_text_sizes(sc_cfg$w)

sc <- ggdraw(sc) +
  draw_label("A", x = sc_cfg$tag_x, y = sc_cfg$y["A"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = sc_cfg$ttl_x, y = sc_cfg$y["A"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["A"] - sc_cfg$sub_off, label = pA_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  ) +
  draw_label("B", x = sc_cfg$tag_x, y = sc_cfg$y["B"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = sc_cfg$ttl_x, y = sc_cfg$y["B"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["B"] - sc_cfg$sub_off, label = pB_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  ) +
  draw_label("C", x = sc_cfg$tag_x, y = sc_cfg$y["C"], size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = sc_cfg$ttl_x, y = sc_cfg$y["C"], size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  annotate("text",
    x = sc_cfg$ttl_x, y = sc_cfg$y["C"] - sc_cfg$sub_off, label = pC_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = txt$subtitle / .pt, colour = "grey30"
  )

ggsave(file.path(RPT, "F01_single_col.pdf"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT, "F01_single_col.png"), sc,
  width = sc_cfg$w, height = sc_cfg$h, units = "mm", dpi = 300
)

pB_comp2 <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp2 <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

dc <- wrap_elements(full = pA + theme(plot.margin = margin(8, 2, 10, 2))) +
  wrap_elements(full = pB_comp2 & theme(plot.margin = margin(4, 2, 2, 2))) +
  wrap_elements(full = pC_comp2 & theme(plot.margin = margin(4, 2, 4, 2))) +
  plot_layout(design = "AB\nAC", widths = c(0.35, 0.65), heights = c(1, 1))

dc <- ggdraw(dc) +
  draw_label("A",
    x = dc_cfg$x_a, y = dc_cfg$y_top,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pA_title,
    x = dc_cfg$x_a + dc_cfg$ttl_off, y = dc_cfg$y_top,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_a + dc_cfg$ttl_off, y = dc_cfg$y_top - dc_cfg$sub_off,
    label = pA_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  ) +
  draw_label("B",
    x = dc_cfg$x_bc, y = dc_cfg$y_top,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pB_title,
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_top,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_top - dc_cfg$sub_off,
    label = pB_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  ) +
  draw_label("C",
    x = dc_cfg$x_bc, y = dc_cfg$y_mid,
    size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label(pC_title,
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_mid,
    size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  annotate("text",
    x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_mid - dc_cfg$sub_off,
    label = pC_subtitle,
    parse = TRUE, hjust = 0, vjust = 1, size = dc_cfg$sub_sz / .pt, colour = "grey30"
  )

ggsave(file.path(RPT, "F01.pdf"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT, "F01.png"), dc,
  width = dc_cfg$w, height = dc_cfg$h, units = "mm", dpi = 300
)

message("F01 main composites done")
