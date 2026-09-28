#!/usr/bin/env Rscript
# S1 Table: folds the data frames the S1a and S1b panels leave in
# c_data/sheets into one workbook, then deletes them. Run after S1a.R and S1b.R.

setwd(here::here())

pacman::p_load(dplyr)

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F00/c_data"
SHEETS <- file.path(DAT, "sheets")
sheet_names <- paste0("panel_", LETTERS[1:14])
sheet_rds <- file.path(SHEETS, paste0(sheet_names, ".rds"))
if (!all(file.exists(sheet_rds))) {
  stop("no panel data in ", SHEETS, "; run S1a.R and S1b.R first")
}
sheets <- setNames(lapply(sheet_rds, readRDS), sheet_names)

bench_plot <- sheets$panel_J |>
  filter(method != "Non_imputed")

build_workbook(
  file.path(DAT, "F00_data.xlsx"),
  title = "S1 Table \u2014 quality control",
  description = "Source data behind S1a Figure (normalization, panels A\u2013G) and S1b Figure (imputation, panels H\u2013N). One sheet per panel.",
  overview_df = data.frame(
    Sheet = sheet_names,
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
  sheet_specs = lapply(names(sheets), function(n) {
    list(name = n, df = as.data.frame(sheets[[n]]))
  })
)
unlink(SHEETS, recursive = TRUE)

message("F00 complete")
