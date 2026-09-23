# 05_supplement_sensitivity — how much proteome variance sits on the
# supplementation stratum, and does adjusting for it move the older adults'
# training response? Companion to 02_supplement_covariate.R, which asks
# whether supplement can enter the design at all; this one measures the
# variance it carries and tests the training delta between arms.
#
# Supplement is read from YvO_meta.xlsx, not YvO_pheno_calc.xlsx. The
# phenotype workbook leaves the arm blank for all fourteen EAA participants;
# the rebuild on 2026-09-08 recovered them into the metadata, and reading the
# old source here would silently drop them from every count below.
#
# Supplement is nested within age group with no overlap, so it is tested
# within the older arm alone.

withr::local_dir(here::here())
pacman::p_load(withr, dplyr, tidyr, stringr, readxl, limma, vegan, readr)

dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
dal_imputed <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")
meta <- read_excel("00_input/YvO_meta.xlsx") |>
  mutate(
    subject_key = str_remove(Col_ID, "_(Pre|Post)$"),
    supplement = factor(supplement)
  ) |>
  mutate(parent_study = if_else(str_detect(subject_key, "^OP_|^YP_"),
    "peanut_protein", "other"
  )) |>
  filter(Col_ID %in% colnames(dal$data))

stopifnot(!anyNA(meta$supplement))

cat("Participants per supplement group by age:\n")
subj <- distinct(meta, subject_key, Group, supplement)
print(with(subj, table(supplement, Group)))
cat(
  "\nSupplement is nested within age group, so the two terms are not",
  "separately estimable in a model carrying both.\n\n"
)

old <- filter(meta, Group == "Old")
old_mat <- dal$data[, old$Col_ID]
old_mat_imputed <- dal_imputed$data[, old$Col_ID]
old$supplement <- droplevels(old$supplement)

cat("Older adults, supplement against parent study:\n")
print(with(
  distinct(old, subject_key, supplement, parent_study),
  table(supplement, parent_study)
))

# Global variance first, with parent study fitted ahead of supplement so the
# supplement term cannot absorb a cohort difference. PERMANOVA needs a
# complete matrix; the non-imputed DAList carries 12.1% missing (limma
# handles that per-protein in the DEP fit below, but adonis2 can't), so this
# step alone uses the imputed matrix, matching how the manuscript's other
# whole-proteome multivariate analyses (PCA, WGCNA) already use it.
set.seed(42)
permanova <- adonis2(
  t(old_mat_imputed) ~ parent_study + supplement + Timepoint,
  data = old, method = "euclidean", permutations = 999, by = "terms"
)
cat("\nPERMANOVA on older adults (Euclidean, 999 permutations):\n")
print(permanova)

# Sequential, so parent_study is fitted first and supplement is measured on
# what remains. Each participant contributes two rows and the permutation is
# unrestricted, which is liberal for the between-subject terms; the result is
# quoted as the share of variance carried, not as a test of the training
# response, which is what the delta model below tests.
write_csv(
  tibble::rownames_to_column(as.data.frame(permanova), "term"),
  "03_DEP/c_data/05_supplement_permanova.csv"
)

# Does the training response differ by supplement, older adults only?
paired <- old |>
  count(subject_key) |>
  filter(n == 2) |>
  pull(subject_key)

delta <- vapply(
  paired,
  function(s) old_mat[, paste0(s, "_Post")] - old_mat[, paste0(s, "_Pre")],
  numeric(nrow(old_mat))
)
supp_of <- droplevels(old$supplement[match(paired, old$subject_key)])

design <- model.matrix(~ 0 + supp_of)
colnames(design) <- levels(supp_of)
fit <- eBayes(lmFit(delta, design))
ftest <- topTable(fit,
  coef = seq_len(ncol(design)), number = Inf,
  sort.by = "F"
)

# Non-imputed delta carries real missingness (a protein missing Pre or Post
# for a subject drops that subject from its own per-protein fit); testable
# below means the F-test was estimable at all for that protein.
testable <- ftest[!is.na(ftest$adj.P.Val), ]
cat(sprintf(
  "\nTraining-delta differences between BRJ / PLA / PP in older adults (n = %d, %d/%d proteins testable):\n",
  length(paired), nrow(testable), nrow(ftest)
))
cat(sprintf("  proteins at FDR < 0.05: %d\n", sum(testable$adj.P.Val < 0.05)))
cat(sprintf("  proteins at FDR < 0.10: %d\n", sum(testable$adj.P.Val < 0.10)))
cat(sprintf(
  "  proteins at p    < 0.05: %d of %d (%.1f%%, chance expectation 5%%)\n",
  sum(testable$P.Value < 0.05), nrow(testable),
  100 * mean(testable$P.Value < 0.05)
))
cat(sprintf("  smallest adjusted p-value: %.3f\n", min(testable$adj.P.Val)))

# BRJ against PLA within the one parent study that randomised both, so the
# comparison is not carrying a cohort difference.
# Arm names come from YvO_meta.xlsx, which spells them out rather than using
# the phenotype workbook's three-letter codes.
brj_pla <- paired[
  supp_of %in% c("Beetroot juice", "Placebo") &
    old$parent_study[match(paired, old$subject_key)] == "other"
]
supp_bp <- droplevels(old$supplement[match(brj_pla, old$subject_key)])
fit_bp <- eBayes(lmFit(delta[, brj_pla], model.matrix(~supp_bp)))
tt_bp <- topTable(fit_bp, coef = 2, number = Inf)
tt_bp_testable <- tt_bp[!is.na(tt_bp$adj.P.Val), ]

cat(sprintf(
  "\nBRJ vs PLA within the beetroot study only (n = %d BRJ, %d PLA, %d/%d proteins testable):\n",
  sum(supp_bp == "Beetroot juice"), sum(supp_bp == "Placebo"), nrow(tt_bp_testable), nrow(tt_bp)
))
cat(sprintf("  proteins at FDR < 0.10: %d\n", sum(tt_bp_testable$adj.P.Val < 0.10)))
cat(sprintf(
  "  proteins at p    < 0.05: %d of %d (%.1f%%)\n",
  sum(tt_bp_testable$P.Value < 0.05), nrow(tt_bp_testable),
  100 * mean(tt_bp_testable$P.Value < 0.05)
))
cat(sprintf("  smallest adjusted p-value: %.3f\n", min(tt_bp_testable$adj.P.Val)))
