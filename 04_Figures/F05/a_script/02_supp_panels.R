setwd(here::here())

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

BASE <- "04_Figures/F05"

message("sourcing F05 supp QC panels")
source("04_Figures/F05/a_script/_supp_qc_soft_threshold.R")
source("04_Figures/F05/a_script/_supp_qc_dendrogram.R")
source("04_Figures/F05/a_script/_supp_qc_compartment.R")
source("04_Figures/F05/a_script/_supp_qc_bicor.R")

message("sourcing F05 supp QC composite")

pacman::p_load(patchwork, cowplot, png, grid)

# ggrepel places labels by a stochastic search, so an unseeded render puts
# them somewhere new each time. run_all.R runs each script in its own
# Rscript child, which starts from a time-seeded RNG.
set.seed(42)

RPT_SRC <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

read_panel <- function(file, dir = RPT_SRC) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel("SUPP_soft_threshold.png")
pB <- read_panel("SUPP_dendrogram.png")
pC <- read_panel("SUPP_compartment_enrichment.png")
pD <- read_panel("SUPP_bicor_sensitivity.png")

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

ggsave(file.path(RPT_PDF, "SUPP_F05_composite.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)
embed_pdf_fonts(file.path(RPT_PDF, "SUPP_F05_composite.pdf"))
ggsave(file.path(RPT_PNG, "SUPP_F05_composite.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm",
  dpi = 300, limitsize = FALSE
)

message("F05 supp QC saved: SUPP_F05_composite.{pdf,png}")

# The triptych, hub-network and preservation panels moved to 01_main_panels.R:
# the workbook is assembled there and three of its sheets come from the CSVs
# they write, so they have to run first.
# Reads the per-module PNGs that 01_main_panels.R renders, so it cannot run any
# earlier than this.
message("sourcing F05 supp module figure")
source("04_Figures/F05/a_script/_supp_mod_figure.R")

message("F05 supp panels complete")
