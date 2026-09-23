#!/usr/bin/env Rscript
# F01 — Phenotype Figure: Master Orchestrator
# One pass: panels -> per-panel stats -> manuscript tables -> composites -> workbook.

withr::local_dir(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

# Panels first: each writes its own source data and a _summary.csv. The table
# script then assembles Table 1B from those summaries rather than refitting,
# so the figure and the table cannot disagree. The supp composite sources the
# two panel scripts it stacks, so they are not listed again here.
source("04_Figures/F01/a_script/05_supp_composite.R")
source("04_Figures/F01/a_script/01_main_panels.R")
source("04_Figures/F01/a_script/04_phenotype_table.R")

# Set after sourcing: each panel script assigns its own DAT, so a value
# defined above would be whatever the last one happened to leave behind.
DAT <- "04_Figures/F01/c_data"

panel_specs <- list(
  list(name = "panel_A_train_volume", path = file.path(DAT, "panel_A_training_volume.csv")),
  list(name = "panel_B_dxa_lbm", path = file.path(DAT, "panel_B_dxa_lbm.csv")),
  list(name = "panel_C_vl_thickness", path = file.path(DAT, "panel_C_vl_thickness.csv")),
  list(name = "supp_deadlift_1rm", path = file.path(DAT, "supp", "panel_A_deadlift_1rm.csv")),
  list(name = "supp_type_II_fcsa", path = file.path(DAT, "supp", "panel_B_type_II_fcsa.csv")),
  list(name = "supp_type_I_fcsa", path = file.path(DAT, "supp", "panel_C_type_I_fcsa.csv")),
  list(name = "supp_dxa_fat_mass", path = file.path(DAT, "supp", "panel_A_dxa_fat_mass.csv")),
  list(name = "supp_fat_to_lean", path = file.path(DAT, "supp", "panel_B_fat_to_lean.csv"))
)

summary_specs <- lapply(
  list.files(DAT, pattern = "_summary\\.csv$", recursive = TRUE),
  \(f) list(name = sub("_summary\\.csv$", "", basename(f)), path = file.path(DAT, f))
)

f01_specs <- c(
  list(
    list(name = "Table_1A_characteristics", df = table_1a),
    list(name = "Table_1B_pre_post", df = table_1b),
    list(name = "Table_1C_composition", df = table_1c),
    list(name = "Per_participant", df = per_subject),
    list(name = "Notes", df = notes)
  ),
  panel_specs
)

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  title = "S2 Table \u2014 participants and training outcomes",
  description = paste(
    "Table 1 as printed, the per-participant values behind it, and the source",
    "data for Figure 1 and S2 Figure."
  ),
  overview_df = data.frame(
    Sheet = c(
      "Table_1A_characteristics", "Table_1B_pre_post", "Table_1C_composition",
      "Per_participant", "Notes",
      "panel_A_train_volume", "panel_B_dxa_lbm", "panel_C_vl_thickness",
      "supp_deadlift_1rm", "supp_type_II_fcsa", "supp_type_I_fcsa",
      "supp_dxa_fat_mass", "supp_fat_to_lean"
    ),
    Description = c(
      "Table 1 upper block: baseline characteristics by age group.",
      "Table 1 lower block: every pre-post measure, with age, time and age-by-time p-values.",
      "How many participants each parent trial and supplement arm contributed, by age group.",
      "One row per participant: every phenotype value behind Table 1.",
      "What each measure is, how it was obtained, and which participants are missing it.",
      "Figure 1A: training volume per participant and group.",
      "Figure 1B: DXA lean body mass pre and post, and the change.",
      "Figure 1C: vastus lateralis thickness pre and post, and the change.",
      "S2 Figure A: deadlift one-repetition maximum.",
      "S2 Figure B: type II fibre cross-sectional area.",
      "S2 Figure C: type I fibre cross-sectional area.",
      "S2 Figure D: whole-body DXA fat mass.",
      "S2 Figure E: fat-to-lean mass ratio."
    ),
    stringsAsFactors = FALSE
  ),
  sheet_specs = f01_specs
)

cleanup_after_workbook(c(panel_specs, summary_specs),
  extra_subdirs = file.path(DAT, "supp")
)

BOX <- Sys.getenv("YVO_BOX_DIR", unset = "")
if (nzchar(BOX) && dir.exists(BOX)) {
  RPT <- "04_Figures/F01/b_reports"
  box_pdf <- file.path(BOX, "02_Figures", "pdf")
  box_png <- file.path(BOX, "02_Figures", "png")
  box_fig_pdf <- file.path(BOX, "03_Supplementary", "figures", "pdf")
  box_fig_png <- file.path(BOX, "03_Supplementary", "figures", "png")
  box_tbl <- file.path(BOX, "03_Supplementary", "tables")
  for (d in c(box_pdf, box_png, box_fig_pdf, box_fig_png, box_tbl)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  file.copy(file.path(RPT, "main/pdf/MAIN_F01_composite.pdf"),
    file.path(box_pdf, "MAIN_F01_composite.pdf"),
    overwrite = TRUE
  )
  file.copy(file.path(RPT, "main/png/MAIN_F01_composite.png"),
    file.path(box_png, "MAIN_F01_composite.png"),
    overwrite = TRUE
  )
  file.copy(file.path(RPT, "supp/pdf/SUPP_F01_phenotypes.pdf"),
    file.path(box_fig_pdf, "S03_Figure_F01.pdf"),
    overwrite = TRUE
  )
  file.copy(file.path(RPT, "supp/png/SUPP_F01_phenotypes.png"),
    file.path(box_fig_png, "S03_Figure_F01.png"),
    overwrite = TRUE
  )
  file.copy(file.path(DAT, "F01_supplementary.xlsx"),
    file.path(box_tbl, "S05_Table_F01.xlsx"),
    overwrite = TRUE
  )
  message("Copied F01 outputs to Box")
}

message("F01 complete")
