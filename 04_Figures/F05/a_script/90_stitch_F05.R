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

message("F05 complete")
