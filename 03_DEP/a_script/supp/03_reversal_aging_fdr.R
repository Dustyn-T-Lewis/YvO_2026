# 03_reversal_aging_fdr: does training in older adults push the age-affected
# proteome back toward young, asked on the FDR-defined Aging set. Runs against
# this pipeline's own combined results and imputed matrix; group-mean axes
# need a complete matrix and the non-imputed DAList carries 12.1% missing.
#
# Aging (Old_Pre - Young_Pre) and Training_Old (Old_Post - Old_Pre) share
# Old_Pre with opposite signs, so their correlation is negative before any
# biology (Oldham 1962; Smyth & Altman 2013), and selecting proteins on Aging
# significance amplifies the artifact through regression to the mean: proteins
# whose Old_Pre was extreme partly by noise are picked, then "revert" by
# construction. Every inferential number here therefore comes from machinery
# that preserves or breaks that structure explicitly.
#
# Selected set is the Aging FDR < 0.05 proteins from this pipeline's own fit;
# the selection inside the permutation and split-half tracks uses the top-n by
# |group-mean aging logFC| so the null can re-select the same way; nulls are
# coupling-preserving within-subject Pre/Post swaps (N = 1000); split-halves
# put the aging axis and the training axis in disjoint older subjects with
# selection computed only in the aging half (N = 500). The limma-coefficient
# correlation on the actual FDR set is reported as descriptive only.

pacman::p_load(withr, readr, dplyr, tibble)

withr::local_dir(here::here())

N_PERM <- 1000
N_SPLIT <- 500

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
dal <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")

meta <- readxl::read_excel("00_input/YvO_meta.xlsx") |>
  dplyr::filter(Col_ID %in% colnames(dal$data))
x <- dal$data[, meta$Col_ID]
subject <- sub("_(Pre|Post)$", "", meta$Col_ID)

old_pre <- meta$Col_ID[meta$Group == "Old" & meta$Timepoint == "Pre"]
old_post <- meta$Col_ID[meta$Group == "Old" & meta$Timepoint == "Post"]
names(old_pre) <- subject[match(old_pre, meta$Col_ID)]
names(old_post) <- subject[match(old_post, meta$Col_ID)]
paired <- intersect(names(old_pre), names(old_post))
old_pre <- old_pre[paired]
old_post <- old_post[paired]
young_pre <- meta$Col_ID[meta$Group == "Young" & meta$Timepoint == "Pre"]

group_mean <- function(cols) rowMeans(x[, cols, drop = FALSE])
young_ref <- group_mean(young_pre)

aging_set <- dep |>
  filter(!is.na(adj.P.Val_Aging), adj.P.Val_Aging < 0.05)
n_sel <- nrow(aging_set)

desc_df <- aging_set |>
  filter(!is.na(logFC_Training_Old))
desc_r <- cor(desc_df$logFC_Aging, desc_df$logFC_Training_Old)
desc_opposite <- mean(sign(desc_df$logFC_Aging) != sign(desc_df$logFC_Training_Old))

sel_top <- function(aging_axis) {
  order(-abs(aging_axis))[seq_len(n_sel)]
}

aging_obs <- group_mean(old_pre) - young_ref
training_obs <- group_mean(old_post) - group_mean(old_pre)
sel_obs <- sel_top(aging_obs)
obs_r <- cor(aging_obs[sel_obs], training_obs[sel_obs])

set.seed(42)
null_r <- replicate(N_PERM, {
  flip <- runif(length(paired)) < 0.5
  pre_s <- ifelse(flip, old_post, old_pre)
  post_s <- ifelse(flip, old_pre, old_post)
  aging_p <- group_mean(pre_s) - young_ref
  training_p <- group_mean(post_s) - group_mean(pre_s)
  sel <- sel_top(aging_p)
  cor(aging_p[sel], training_p[sel])
})
p_coupling <- (1 + sum(null_r <= obs_r)) / (1 + N_PERM)

set.seed(42)
half <- length(paired) %/% 2
split_r <- replicate(N_SPLIT, {
  o <- sample(length(paired))
  a <- o[seq_len(half)]
  b <- o[(half + 1):length(paired)]
  aging_a <- group_mean(old_pre[a]) - young_ref
  training_b <- group_mean(old_post[b]) - group_mean(old_pre[b])
  sel <- sel_top(aging_a)
  cor(aging_a[sel], training_b[sel])
})

reversal_aging <- tibble(
  estimator = c(
    "limma coefficients on the FDR set (descriptive)",
    "Group-mean axes, top-n selection, observed",
    "Coupling-preserving null with re-selection",
    "Split-half: selection and aging in half A, training in half B"
  ),
  controls_coupling = c(FALSE, FALSE, TRUE, TRUE),
  controls_selection = c(FALSE, FALSE, TRUE, TRUE),
  r = c(desc_r, obs_r, mean(null_r), mean(split_r)),
  lo = c(NA, NA, quantile(null_r, 0.025), quantile(split_r, 0.025)),
  hi = c(NA, NA, quantile(null_r, 0.975), quantile(split_r, 0.975)),
  p_value = c(NA, NA, p_coupling, NA),
  n_proteins = n_sel,
  n_subjects = c(NA, length(paired), length(paired), length(paired)),
  pct_opposite_sign = c(round(100 * desc_opposite, 1), NA, NA, NA)
) |>
  mutate(across(c(r, lo, hi), \(v) round(v, 3)))

write_csv(reversal_aging, "03_DEP/c_data/03_reversal_aging_fdr.csv")

# The draws themselves, not just their summary: the figure that answers the
# circularity objection shows the observed value sitting inside this null, and
# a reader cannot see that from an interval.
write_csv(
  tibble(estimator = rep(c("coupling_null", "split_half"),
                         c(length(null_r), length(split_r))),
         r = c(null_r, split_r)),
  "03_DEP/c_data/03_reversal_null_draws.csv"
)

print(as.data.frame(reversal_aging))
