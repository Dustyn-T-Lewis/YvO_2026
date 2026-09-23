#!/usr/bin/env Rscript
# F06 — Phenotype Prediction: Master Orchestrator
# Sources main panels (all data-generating scripts + composite + xlsx),
# then supp panels (reads pre-rendered PNGs, builds composites).

setwd(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F06/c_data"

message("F06: Running main panels + composite + xlsx")
source("04_Figures/F06/a_script/01_main_panels.R")

message("F06: Running supp composite")
source("04_Figures/F06/a_script/02_supp_panels.R")

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT, pattern = "\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

# Copy to Box manuscript directory
BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
if (nzchar(BOX) && dir.exists(BOX)) {
  RPT <- "04_Figures/F06/b_reports"
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F06_composite.pdf"),
            file.path(box_pdf, "MAIN_F06_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F06_composite.png"),
            file.path(box_png, "MAIN_F06_composite.png"), overwrite = TRUE)
  # S7 Figure (the supplementary composite)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F06_composite.pdf"),
            file.path(box_fig_pdf, "S09_Figure_F06.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F06_composite.png"),
            file.path(box_fig_png, "S09_Figure_F06.png"), overwrite = TRUE)
  # S7 Table
  file.copy(file.path(DAT, "F06_supplementary.xlsx"),
            file.path(box_tbl, "S11_Table_F06.xlsx"), overwrite = TRUE)
  message("Copied F06 outputs to Box")
}

message("F06 complete")
