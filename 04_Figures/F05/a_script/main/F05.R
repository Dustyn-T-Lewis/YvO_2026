#!/usr/bin/env Rscript
# Figure 5: the module-trait heatmap and the module NES scatters, placed as
# rasters on one canvas and cropped.

setwd(here::here())

source("04_Figures/shared/style.R")

pacman::p_load(patchwork, cowplot, png, grid, dplyr, tidyr)

BASE <- "04_Figures/F05"

PANELS <- "04_Figures/F05/a_script/main/panels"
source_panel(file.path(PANELS, "A_module_trait_heatmap.R"))
source_panel(file.path(PANELS, "B_nes_scatters.R"))

RPT <- file.path(BASE, "b_reports", "main")

PANEL_DIR <- file.path(BASE, "b_reports", "main", "panels")
read_panel <- function(file, dir = PANEL_DIR) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob <- read_panel("A_module_trait_heatmap.png")
pB_grob <- read_panel("B_nes_scatters.png")
pB_leg_grob <- read_panel("B_nes_scatters_legend.png")

COMP_W <- 470
COMP_H <- 300
# Panel B's box is solved, not chosen, so the scatters' outer panel borders
# land on the heatmap's top and bottom tile borders. Both panels are rasters
# fitted to their box with the aspect preserved, so the heatmap borders sit at
# canvas y 0.751686 and 0.258944, and the scatter borders at 0.056307 and
# 0.990219 of panel B's own height. Re-measure both if either panel's internal
# spacing changes: panel B's axis type sets where its borders fall.
GRID_TOP <- 0.756847
GRID_BOT <- 0.229236
CROP_L <- 0.01
CROP_R <- 0.80
CROP_B <- 0.17
CROP_T <- 0.85
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

# Type scales against SAVE_W, not COMP_W: the composite is drawn on a 470 mm
# canvas but cropped to about 371 mm before saving, and it is the saved width
# that shrinks to the print width. Sizing off COMP_W overshoots by ~27% and
# pushes panel B's title past the right edge.
f06_txt <- composite_text_sizes(SAVE_W, 178)

# Tag, title and subtitle for each panel, positioned in canvas mm.
header <- function(tag, title, sub, let_x, ttl_x) {
  list(
    draw_label(tag,
      x = mm2x(let_x), y = mm2y(248),
      size = f06_txt$tag, fontface = "bold", hjust = 0, vjust = 1
    ),
    draw_label(title,
      x = mm2x(ttl_x), y = mm2y(248),
      size = f06_txt$title, fontface = "bold", hjust = 0, vjust = 1
    ),
    draw_label(sub,
      x = mm2x(ttl_x), y = mm2y(243),
      size = f06_txt$subtitle, fontface = "bold.italic", colour = "grey40",
      hjust = 0, vjust = 1
    )
  )
}

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B first (behind) so its white left margin is hidden under Panel A
  draw_grob(pB_grob,
    x = 0.43, y = GRID_BOT, width = 0.55, height = GRID_TOP - GRID_BOT,
    hjust = 0, vjust = 0
  ) +
  # Panel B's legend is clipped: its three columns are wider than the gap
  # between the heatmap and the crop edge, so the first column loses a few
  # characters. Widening it only moves the clipping to the right-hand column;
  # the fix is to re-flow the legend, which is layout surgery.
  draw_grob(pB_leg_grob,
    x = 0.70, y = 0.195, width = 0.24, height = 0.04,
    hjust = 0.5, vjust = 0.5
  ) +
  draw_grob(pA_grob,
    x = 0.01, y = 0, width = 0.60, height = 0.96,
    hjust = 0, vjust = 0
  ) +
  header(
    "A", "WGCNA Module-Trait Associations",
    paste0(
      "9 modules, labelled enrichment then hub | ",
      "LMM (BH per-trait) | * BH, dashed nominal, bold |r| >= 0.40"
    ),
    let_x = 15, ttl_x = 58
  ) +
  header(
    "B", "Module-Level NES Scatters", "fGSEA on module-member t-stats",
    let_x = 285, ttl_x = 294
  )

ggsave(file.path(RPT, "F05.pdf"), composite_final,
  width = SAVE_W, height = SAVE_H, units = "mm",
  device = get_raster_pdf_device()
)
embed_pdf_fonts(file.path(RPT, "F05.pdf"))
ggsave(file.path(RPT, "F05.png"), composite_final,
  width = SAVE_W, height = SAVE_H, units = "mm",
  dpi = 300
)
message("F05 main composite done")
