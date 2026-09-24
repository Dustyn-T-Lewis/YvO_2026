#!/usr/bin/env Rscript
# S2 Figure: the five phenotype pre/post pairs in one figure.
#
# Strength/CSA and body composition used to render as two figures on two
# canvases, 88 mm over three rows against 110 mm over two, so the same panel
# was drawn half again as tall depending on which file it landed in. One
# canvas at the tighter row height ends that.

setwd(here::here())

pacman::p_load(withr, patchwork, cowplot)

source("04_Figures/shared/style.R")

# Tag to panel script, in reading order. A and D set the seed; B, C and E
# draw their jitter from the stream the panel before them leaves.
PANELS <- "04_Figures/F01/a_script/supp/panels"
ROW_SCRIPT <- c(
  A = "S2_A_deadlift_1rm.R", B = "S2_B_type_II_fcsa.R",
  C = "S2_C_type_I_fcsa.R", D = "S2_D_dxa_fat_mass.R", E = "S2_E_fat_to_lean.R"
)
rows <- lapply(ROW_SCRIPT, \(f) {
  p <- source_panel(file.path(PANELS, f))
  list(
    left = strip_for_composite(p[[1]]), right = strip_for_composite(p[[2]]),
    title = p[[1]]$labels$title, subtitle = attr(p, "subtitle_expr")
  )
})

RPT <- "04_Figures/F01/b_reports/supp"

N_ROWS <- length(rows)
COMP_W <- 85
ROW_H <- 32
COMP_H <- ROW_H * N_ROWS
PT_PER_MM <- 72 / 25.4
# The top margin doubles as the band the heading is drawn into, so it has to
# hold a 5 pt title over a 5 pt subtitle with a little air.
HEAD_PT <- 15
SIDE_PT <- 2

composite <- wrap_plots(
  lapply(rows, \(r) list(r$left, r$right)) |> unlist(recursive = FALSE),
  ncol = 2, widths = c(0.65, 0.35),
  heights = rep(1, N_ROWS)
) &
  theme(plot.margin = margin(HEAD_PT, SIDE_PT, SIDE_PT, SIDE_PT))

txt <- composite_text_sizes(COMP_W)
TAG_SZ <- txt$tag - 3
TTL_SZ <- txt$title - 2
SUB_SZ <- txt$subtitle - 1
X_TAG <- 0.056
X_TTL <- 0.116

# `&` hands the margin to the assembly as well as to each panel, so the five
# rows share what is left of the canvas after the outer margin rather than
# taking a clean fifth each. Placing headings at i/N without that correction
# is what walked them onto the panel frame a little further down every row.
CANVAS_PT <- COMP_H * PT_PER_MM
ROW_PT <- (CANVAS_PT - HEAD_PT - SIDE_PT) / N_ROWS
TTL_OFF_PT <- 1.5
# 5.7 pt of leading under a 5 pt title: any tighter and the two lines share
# bounding boxes, which is what the ascender of the subtitle was doing.
SUB_OFF_PT <- 7.2
SUB_OFF <- (SUB_OFF_PT - TTL_OFF_PT) / CANVAS_PT

composite <- ggdraw(composite)
for (i in seq_along(rows)) {
  r <- rows[[i]]
  row_top <- HEAD_PT + (i - 1) * ROW_PT
  y <- 1 - (row_top + TTL_OFF_PT) / CANVAS_PT
  composite <- composite +
    draw_label(names(rows)[i],
      x = X_TAG, y = y, size = TAG_SZ,
      fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(r$title,
      x = X_TTL, y = y, size = TTL_SZ,
      fontface = "bold", hjust = 0, vjust = 1
    ) +
    # annotate() rather than draw_label(): the subtitle is plotmath, so only
    # the significant RM-ANOVA terms come out bold, and draw_label cannot
    # parse an expression.
    annotate("text",
      x = X_TTL, y = y - SUB_OFF, label = r$subtitle, parse = TRUE,
      hjust = 0, vjust = 1, size = SUB_SZ / .pt, colour = "grey30"
    )
}

ggsave(file.path(RPT, "S2.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT, "S2.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

caption_supp(composite, "S2", COMP_W, COMP_H, RPT)

message("S2 done")
