#!/usr/bin/env Rscript
# S6 Figure: network construction diagnostics. Soft threshold, dendrogram,
# compartment enrichment and the bicor comparison, placed as rasters.

setwd(here::here())

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

BASE <- "04_Figures/F05"

PANELS <- "04_Figures/F05/a_script/supp/panels"
for (f in c("S6_A_soft_threshold.R", "S6_B_dendrogram.R", "S6_C_compartment.R",
            "S6_D_bicor.R")) {
  source_panel(file.path(PANELS, f))
}

pacman::p_load(patchwork, cowplot, png, grid)

# ggrepel places labels by a stochastic search, so an unseeded render puts
# them somewhere new each time. run_all.R runs each script in its own
# Rscript child, which starts from a time-seeded RNG.
set.seed(42)

RPT_SRC <- file.path(BASE, "b_reports", "supp", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp")
RPT_PNG <- file.path(BASE, "b_reports", "supp")

read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel("S6_A_soft_threshold.png")
pB <- read_panel("S6_B_dendrogram.png")
pC <- read_panel("S6_C_compartment.png")
pD <- read_panel("S6_D_bicor.png")

bottom_row <- wrap_elements(full = pC) + wrap_elements(full = pD) +
  plot_layout(widths = c(1, 1))

composite <- wrap_elements(full = pA) /
  wrap_elements(full = pB) /
  bottom_row +
  plot_layout(heights = c(0.30, 0.35, 0.35)) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(4, 6, 4, 6)
    )
  )

COMP_W <- 250
COMP_H <- 330

TAG_SZ <- composite_text_sizes(COMP_W, 178)$tag

composite <- ggdraw(composite) +
  draw_label("A", x = 0.02, y = 0.960, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.02, y = 0.660, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.02, y = 0.320, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.52, y = 0.320, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

pdf_device <- get_raster_pdf_device()

ggsave(file.path(RPT_PDF, "S6.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)
embed_pdf_fonts(file.path(RPT_PDF, "S6.pdf"))
ggsave(file.path(RPT_PNG, "S6.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm",
  dpi = 300, limitsize = FALSE
)

caption_supp(composite, "S6", COMP_W, COMP_H, RPT_PDF)

message("S6 done")
