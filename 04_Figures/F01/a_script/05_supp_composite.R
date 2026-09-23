#!/usr/bin/env Rscript
# F01 Supp composite — the five phenotype pre/post pairs in one figure.
#
# Strength/CSA and body composition used to render as two figures on two
# canvases, 88 mm over three rows against 110 mm over two, so the same panel
# was drawn half again as tall depending on which file it landed in. One
# canvas at the tighter row height ends that.

withr::local_dir(here::here())

pacman::p_load(withr, patchwork, cowplot)

source("04_Figures/shared/style.R")

source("04_Figures/F01/a_script/02_supp_panels.R")
source("04_Figures/F01/a_script/03_body_comp_panels.R")

RPT_PNG <- "04_Figures/F01/b_reports/supp/png"
RPT_PDF <- "04_Figures/F01/b_reports/supp/pdf"

# Tag to panel prefix, in reading order. The two panel scripts publish their
# plots under these prefixes, so get() here is the other half of the assign()
# the template ends on.
ROW_PREFIX <- c(A = "pSA", B = "pSB", C = "pSC", D = "pFA", E = "pFB")
rows <- lapply(ROW_PREFIX, \(pfx) {
  parts <- c(
    left = "_left", right = "_right",
    title = "_title", subtitle = "_subtitle"
  )
  lapply(parts, \(suffix) get(paste0(pfx, suffix)))
})

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

ggsave(file.path(RPT_PDF, "SUPP_F01_phenotypes.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device()
)
ggsave(file.path(RPT_PNG, "SUPP_F01_phenotypes.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message("F01 supp composite done")
