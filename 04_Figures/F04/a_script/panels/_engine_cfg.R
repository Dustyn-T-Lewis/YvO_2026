# The loaded objects the trajectory and fry engines take, shared by panels C
# and F.

pacman::p_load(readr, readxl)

# Shared engines take loaded objects, not paths, so the limma tree can supply
# the same structures from its own DAList without the engines knowing which
# pipeline they are running under.
# limma drops outliers upstream, so the metadata sheet carries samples the
# imputed matrix no longer has. Subset to the intersection before anything
# downstream builds a design from it.
abundance <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")$data
meta <- as.data.frame(read_excel("00_input/YvO_meta.xlsx"))
meta <- meta[meta$Col_ID %in% colnames(abundance), ]

engine_cfg <- list(
  meta = meta,
  matrix = as.matrix(abundance[, meta$Col_ID]),
  dep_df = read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE),
  group_time = factor(meta$Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
  )
)
