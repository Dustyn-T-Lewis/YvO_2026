#!/usr/bin/env Rscript
# Copies the manuscript files into Supplementary/ under the names the journal
# and the manuscript use: Figure_1.pdf, S1a_Figure.pdf (with its legend),
# S1_Table.xlsx. Runs last; every file it copies is written by an earlier step.

setwd(here::here())

OUT <- "Supplementary"

figures <- setNames(
  sprintf("04_Figures/F0%d/b_reports/main/F0%d.pdf", 1:6, 1:6),
  sprintf("Figure_%d.pdf", 1:6)
)

supp_figures <- list.files(
  c("03_DEP/b_reports/supp", sprintf("04_Figures/F0%d/b_reports/supp", 0:6)),
  pattern = "^S[0-9]+[ab]?_Figure\\.pdf$", full.names = TRUE
)
supp_figures <- setNames(supp_figures, basename(supp_figures))

supp_tables <- setNames(
  c(
    sprintf("04_Figures/F0%d/c_data/F0%d_data.xlsx", 0:6, 0:6),
    "01_normalization/c_data/01_normalization.xlsx",
    "02_imputation/c_data/02_imputation.xlsx",
    "03_DEP/c_data/03_DEP_results.xlsx"
  ),
  sprintf("S%d_Table.xlsx", 1:10)
)

sets <- list(
  "Figures_1-6" = figures,
  "Supplementary_Figures" = supp_figures,
  "Supplementary_Tables" = supp_tables
)
stopifnot(
  "expected 12 supplementary figures" = length(supp_figures) == 12,
  "missing source files; run the pipeline first" = file.exists(unlist(sets))
)

unlink(OUT, recursive = TRUE)
for (dir in names(sets)) {
  dir.create(file.path(OUT, dir), recursive = TRUE)
  file.copy(sets[[dir]], file.path(OUT, dir, names(sets[[dir]])))
}
message("Supplementary/: ", paste(lengths(sets), names(sets), collapse = ", "))
