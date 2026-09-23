# F05 master orchestrator: run YvO_WGCNA_run.R separately first.


setwd(here::here())

source("04_Figures/F05/a_script/01_main_panels.R")
source("04_Figures/F05/a_script/02_supp_panels.R")

# After both, so the four QC sheets hold this run's data rather than the last.
build_workbook(
  file.path("04_Figures/F05", "c_data", "F05_supplementary.xlsx"),
  title = "S6 Table \u2014 co-expression network",
  description = "Source data for Figure 5, S6 Figure and S8 Figure: module assignments, eigengenes, module\u2013trait associations, the evidence behind each module name, and network construction diagnostics.",
  overview_df = f06_overview,
  sheet_specs = f06_specs
)
f06_cleanup()

BASE <- "04_Figures/F05"
BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
RPT <- file.path(BASE, "b_reports")

if (nzchar(BOX) && dir.exists(BOX)) {
  box_pdf     <- file.path(BOX, "02_Figures", "pdf")
  box_png     <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl     <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(RPT, "main/pdf/MAIN_F05_composite.pdf"),
            file.path(box_pdf, "MAIN_F05_composite.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "main/png/MAIN_F05_composite.png"),
            file.path(box_png, "MAIN_F05_composite.png"), overwrite = TRUE)
  # S7 Figure (network diagnostics + module triptychs)
  file.copy(file.path(RPT, "supp/pdf/SUPP_F05_composite.pdf"),
            file.path(box_fig_pdf, "S08_Figure_F05.pdf"), overwrite = TRUE)
  file.copy(file.path(RPT, "supp/png/SUPP_F05_composite.png"),
            file.path(box_fig_png, "S08_Figure_F05.png"), overwrite = TRUE)
  # S7 Figure continued: module triptychs (individual pages)
  triptych_pdfs <- list.files(file.path(RPT, "supp/pdf/modules"),
                              pattern = "^SUPP_triptych_", full.names = TRUE)
  for (tp in triptych_pdfs)
    file.copy(tp, file.path(box_fig_pdf, basename(tp)), overwrite = TRUE)
  # S9 Table
  file.copy(file.path(BASE, "c_data/F05_supplementary.xlsx"),
            file.path(box_tbl, "S10_Table_F05.xlsx"), overwrite = TRUE)
  message("Copied F05 outputs to Box")
}

message("F05 complete")
