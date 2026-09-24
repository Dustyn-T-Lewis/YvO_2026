#!/usr/bin/env Rscript
# S4b Figure: the four significance-criterion heatmaps, one per page.

setwd(here::here())

pacman::p_load(ggplot2, qpdf)

source("04_Figures/shared/style.R")

PANELS <- "04_Figures/F03/a_script/supp/panels"
RPT <- "04_Figures/F03/b_reports/supp"
pages <- c(A = "S4b_A_aging_fdr.R", B = "S4b_B_aging_pi.R",
           C = "S4b_C_training_fdr.R", D = "S4b_D_training_pi.R")

page_pdfs <- character()
for (tag in names(pages)) {
  p <- source_panel(file.path(PANELS, pages[[tag]]))
  h <- attr(p, "height_mm")
  ggsave(file.path(RPT, sprintf("S4b_%s.png", tag)), p,
    width = 178, height = h, units = "mm", dpi = 300, bg = "white",
    device = ragg::agg_png
  )
  page_pdfs[tag] <- tempfile(fileext = ".pdf")
  ggsave(page_pdfs[tag], p,
    width = 178, height = h, units = "mm",
    device = get_pdf_device(), bg = "white"
  )
}
# The heatmaps grow with their protein lists, so each page has its own height
# and one device cannot draw them all.
pdf_combine(page_pdfs, file.path(RPT, "S4b.pdf"))

caption_supp(file.path(RPT, "S4b.pdf"), "S4b", 178, NA, RPT)

message("S4b done")
