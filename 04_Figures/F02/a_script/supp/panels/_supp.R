# Setup shared by the S3 Figure panels: packages, the seed, the shared theme,
# the F02 overrides and the normalized data.

pacman::p_load(dplyr, tidyr, stringr, readr, readxl, ggplot2, ggrepel, ggbeeswarm, patchwork, cowplot)

# ggrepel places labels by a stochastic search, so an unseeded render puts
# them somewhere new each time. run_all.R runs each script in its own
# Rscript child, which starts from a time-seeded RNG.
set.seed(42)

source("04_Figures/shared/style.R")

# F02-specific overrides
HEATMAP_LO <- "#2166AC"; HEATMAP_HI <- "#B2182B"
BASE_COUNT <- BASE_COUNT + 1.0
BASE_GENE  <- BASE_GENE  + 0.8
BASE_STAT  <- BASE_STAT  + 0.5

RPT_PNG <- "04_Figures/F02/b_reports/supp/panels"
RPT_PDF <- RPT_PNG
DAT     <- "04_Figures/F02/c_data"
for (d in c(RPT_PNG, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()

dal_norm <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
norm_mat <- as.matrix(dal_norm$data)
norm_meta <- as_tibble(dal_norm$metadata) |>
  mutate(age     = factor(Group, levels = c("Young", "Old")),
         time    = factor(Timepoint, levels = c("Pre", "Post")),
         subject = sub("_(Pre|Post)$", "", Col_ID),
         group   = factor(Group_Time,
           levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post"))) |>
  rename(sample_id = Col_ID)

samp_names <- norm_meta$sample_id

# Annotation from DAList
ann_df <- as_tibble(dal_norm$annotation) |>
  select(uniprot_id, gene, protein, description)
norm_df <- bind_cols(ann_df, as_tibble(norm_mat))
