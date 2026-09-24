# Setup shared by the Figure 2 panels: packages, the shared theme, the F02
# overrides and the inputs panels A to C read.

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, readxl, ggplot2, patchwork, cowplot, vegan, ComplexHeatmap, purrr)

source("04_Figures/shared/style.R")

# F02-specific overrides (from the old F02/style.R)
HEATMAP_LO <- "#2166AC"
HEATMAP_HI <- "#B2182B"
BASE_COUNT <- BASE_COUNT + 1.0
BASE_GENE <- BASE_GENE + 0.8
BASE_STAT <- BASE_STAT + 0.5

# One size for every in-panel stat box (A's PERMANOVA label, B's per-contrast
# labels) so the two cannot drift apart when either panel is retuned.
STAT_BOX_PT <- 4.4

PNL_PNG <- "04_Figures/F02/b_reports/main/panels"
PNL_PDF <- PNL_PNG
DAT <- "04_Figures/F02/c_data"
for (d in c(PNL_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
  show_col_types = FALSE
)

dal_imp <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")
imp_mat <- as.matrix(dal_imp$data)
imp_meta <- as_tibble(dal_imp$metadata) |>
  mutate(
    age = factor(Group, levels = c("Young", "Old")),
    time = factor(Timepoint, levels = c("Pre", "Post")),
    group = factor(Group_Time,
      levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
    ),
    subject = sub("_(Pre|Post)$", "", Col_ID)
  ) |>
  rename(sample_id = Col_ID)

DEP_XLSX <- "03_DEP/c_data/03_DEP_results.xlsx"
CONTRASTS <- c("Aging", "Training_Young", "Training_Old", "Interaction")

SET_LABELS <- c(
  Aging = "Aging", Training_Young = "Tr.(Y)",
  Training_Old = "Tr.(O)", Interaction = "Inter."
)
