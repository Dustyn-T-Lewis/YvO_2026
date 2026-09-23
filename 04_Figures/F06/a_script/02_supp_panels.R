#!/usr/bin/env Rscript
# F06 Supplementary Composite Stitch
#
# Reads pre-rendered PNGs from b_reports/supp/png/panels/ and composites them.
#
# Output 1: SUPP_F06_composite.pdf  — panels A and B above, C centred below
# Output 2: SUPP_F06_composite.png  — the same composite
#
# Panel sources (all in b_reports/supp/png/panels/):
#   A = Per-module ROC grid      (SUPP_F06_module_grid.png)
#   B = Panel B full sweep       (SUPP_F06_panel_B_grid.png)
#
# Three panels were dropped: the two leave-one-subject-out sensitivity plots
# and the multivariate classifier decomposition. Each plotted a handful of
# AUCs that ship as loso_auc_summary, loso_wgcna_refit_summary and
# panel_A_classifier_auc in F06_supplementary.xlsx, where they read as a
# table. Their scripts still run, because those sheets are what they write.

setwd(here::here())

source("04_Figures/shared/style.R")

pacman::p_load(patchwork, cowplot, png, grid)

BASE       <- "04_Figures/F06"
RPT_PNG    <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF    <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PANELS <- file.path(RPT_PNG, "panels")

read_panel <- function(file) {
  path <- file.path(RPT_PANELS, file)
  if (!file.exists(path)) stop("Missing: ", path)
  img <- readPNG(path)
  list(grob = rasterGrob(img, interpolate = TRUE),
       aspect = dim(img)[2] / dim(img)[1])   # width/height
}

# Page 1 panels (A + B, aspect-matched)
pA <- read_panel("SUPP_F06_module_grid.png")
pB <- read_panel("SUPP_F06_panel_B_grid.png")

COMP_H <- 166   # mm tall
wA_mm  <- COMP_H * pA$aspect
wB_mm  <- COMP_H * pB$aspect
COMP_W <- wA_mm + wB_mm

TAG_SZ <- composite_text_sizes(COMP_W, 178)$tag

page1 <- (wrap_elements(full = pA$grob) |
           wrap_elements(full = pB$grob)) +
  plot_layout(widths = c(wA_mm, wB_mm)) &
  theme(plot.margin = margin(0, 0, 0, 0))

TAG_X_A <- 5 / COMP_W
TAG_X_B <- (wA_mm + 5) / COMP_W
TAG_Y   <- 0.992

page1_final <- ggdraw(page1 & theme(plot.margin = margin(2, 2, 2, 2))) +
  draw_label("A", x = TAG_X_A, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = TAG_X_B, y = TAG_Y, size = TAG_SZ,
             fontface = "bold", hjust = 0, vjust = 1)

# Write outputs
graphics.off()
pdf_device <- get_raster_pdf_device()

ggsave(file.path(RPT_PDF, "SUPP_F06_composite.pdf"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
embed_pdf_fonts(file.path(RPT_PDF, "SUPP_F06_composite.pdf"))

ggsave(file.path(RPT_PNG, "SUPP_F06_composite.png"), page1_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message(sprintf("Wrote SUPP_F06_composite (%.0f x %.0f mm)", COMP_W, COMP_H))
