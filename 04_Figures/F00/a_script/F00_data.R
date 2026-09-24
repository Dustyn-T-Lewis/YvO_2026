#!/usr/bin/env Rscript
# S1 Table: folds the data frames the S1a and S1b panels leave in
# c_data/sheets into one workbook, then deletes them. Run after S1a.R and S1b.R.

setwd(here::here())

pacman::p_load(dplyr)

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F00/c_data"
SHEETS <- file.path(DAT, "sheets")
sheet_rds <- file.path(SHEETS, paste0("panel_", LETTERS[1:14], ".rds"))
if (!all(file.exists(sheet_rds))) {
  stop("no panel data in ", SHEETS, "; run S1a.R and S1b.R first")
}
sheets <- setNames(lapply(sheet_rds, readRDS), paste0("panel_", LETTERS[1:14]))

bench_plot <- sheets$panel_J |>
  filter(method != "Non_imputed")

build_workbook(
  file.path(DAT, "F00_data.xlsx"),
  title = "S1 Table \u2014 quality control",
  description = "Source data behind S1a Figure (normalization, panels A\u2013G) and S1b Figure (imputation, panels H\u2013N). One sheet per panel.",
  overview_df = data.frame(
    Sheet = paste0("panel_", LETTERS[1:14]),
    Description = c(
      "Filter cascade",             "Protein missingness histogram",
      "Pre-norm PCA scores",        "Post-norm PCA scores",
      "Eta-squared values",         "Outlier consensus diagnostics",
      "Per-sample missingness",     "Missingness classification",
      "Classification counts",      sprintf("Benchmark ranking (%d methods)", nrow(bench_plot)),
      "Imputation density summary", "MNAR shift audit",
      "Sample integrity",           "DEP counts per contrast"
    )
  ),
  sheet_specs = list(
    list(name = "panel_A", df = as.data.frame(sheets$panel_A)),
    list(name = "panel_B", df = as.data.frame(sheets$panel_B)),
    list(name = "panel_C", df = as.data.frame(sheets$panel_C)),
    list(name = "panel_D", df = as.data.frame(sheets$panel_D)),
    list(name = "panel_E", df = as.data.frame(sheets$panel_E)),
    list(name = "panel_F", df = as.data.frame(sheets$panel_F)),
    list(name = "panel_G", df = as.data.frame(sheets$panel_G)),
    list(name = "panel_H", df = as.data.frame(sheets$panel_H)),
    list(name = "panel_I", df = as.data.frame(sheets$panel_I)),
    list(name = "panel_J", df = as.data.frame(sheets$panel_J)),
    list(name = "panel_K", df = as.data.frame(sheets$panel_K)),
    list(name = "panel_L", df = as.data.frame(sheets$panel_L)),
    list(name = "panel_M", df = as.data.frame(sheets$panel_M)),
    list(name = "panel_N", df = as.data.frame(sheets$panel_N))
  )
)
unlink(SHEETS, recursive = TRUE)

message("F00 complete")
