#!/usr/bin/env Rscript
# F06 Supplementary — Age Discrimination via Multi-Classifier Comparison
# Writes: c_data/panel_A_{multi_classifier_auc,feature_stability,permutation,roc_curves}.csv
#
# 4 classifiers compared head-to-head:
#   1. Pre-only:  baseline eigengenes
#   2. Post-only: post-training eigengenes
#   3. delta-ME:  training response (Post - Pre)
#   4. Combined:  mean(Pre, Post) per subject — AVERAGED not stacked
#
# Sourced by F06_data.R — expects style.R + figure_supplement_helpers.R
# already loaded.

pacman::p_load(tidyverse, pROC)
source("04_Figures/F06/a_script/_loocv.R")

BASE     <- "04_Figures/F06"
DAT      <- "04_Figures/F05/c_data"
DAT_OUT  <- file.path(BASE, "c_data")
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

F05_SUPP <- "04_Figures/F05/c_data/F05_data.xlsx"
stopifnot("run 04_Figures/F05/a_script/F05_data.R first: missing F05_data.xlsx" =
  file.exists(F05_SUPP))

MEs       <- read_matrix_sheet(F05_SUPP, "MEs",     "sample_id")
me_pre    <- read_matrix_sheet(F05_SUPP, "me_pre",  "subject_key")
me_post   <- read_matrix_sheet(F05_SUPP, "me_post", "subject_key")
subj_age  <- read_sheet_df(F05_SUPP, "metadata_subj_age")
common_subj <- read_vector_sheet(F05_SUPP, "common_subj")

pdf_device <- get_pdf_device()

message("SUPP multivariate: multi-classifier age discrimination...")

me_delta <- as.matrix(me_post[common_subj, ]) - as.matrix(me_pre[common_subj, ])

true_labels <- ifelse(
  subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0
)
n_subj  <- length(common_subj)
n_young <- sum(true_labels == 0)
n_old   <- sum(true_labels == 1)

# Run all 4 classifiers
set.seed(42)

clf_names <- c("Pre", "Post", "DeltaME", "Combined")
clf_colors <- c(Pre = "#D6604D", Post = "#5DA5DA",
                DeltaME = "#4CAF50", Combined = "#9B7FBF")

message("  Running Pre-only classifier...")
res_pre  <- run_topk_loocv(true_labels, me_pre[common_subj, ])

message("  Running Post-only classifier...")
res_post <- run_topk_loocv(true_labels, me_post[common_subj, ])

message("  Running delta-ME classifier...")
res_delta <- run_topk_loocv(true_labels, me_delta)

message("  Running Combined (averaged Pre+Post) classifier...")
me_avg <- (as.matrix(me_pre[common_subj, ]) + as.matrix(me_post[common_subj, ])) / 2
res_combined <- run_topk_loocv(true_labels, me_avg)

all_res <- list(Pre = res_pre, Post = res_post,
                DeltaME = res_delta, Combined = res_combined)

clf_summary <- map_dfr(clf_names, function(nm) {
  probs <- all_res[[nm]]$probs
  roc_obj <- tryCatch(roc(true_labels, probs, quiet = TRUE),
                       error = function(e) NULL)
  if (is.null(roc_obj)) {
    return(tibble(classifier = nm, auc = 0.5,
                  ci_lo = NA_real_, ci_hi = NA_real_))
  }
  ci_obj <- ci.auc(roc_obj)
  acc <- mean((probs > 0.5) == true_labels)
  tibble(classifier = nm, auc = as.numeric(auc(roc_obj)),
         ci_lo = ci_obj[1], ci_hi = ci_obj[3],
         accuracy = acc, n = n_subj, n_young = n_young, n_old = n_old)
})

best_clf <- clf_summary |> arrange(desc(auc)) |> slice(1) |> pull(classifier)
best_auc <- clf_summary |> filter(classifier == best_clf) |> pull(auc)

message(sprintf("  AUC summary:"))
for (i in seq_len(nrow(clf_summary))) {
  row <- clf_summary[i, ]
  star <- if (row$classifier == best_clf) " *" else ""
  message(sprintf("    %s: AUC=%.3f [%.3f-%.3f] Acc=%.1f%%%s",
                  row$classifier, row$auc, row$ci_lo, row$ci_hi,
                  row$accuracy * 100, star))
}

feat_stab <- map_dfr(clf_names, function(nm) {
  freq <- c(table(factor(unlist(all_res[[nm]]$selected), levels = colnames(me_pre))))
  tibble(classifier = nm, feature = names(freq),
         n_selected = as.integer(freq),
         pct_selected = round(100 * freq / n_subj, 1))
})
write_csv(feat_stab, file.path(DAT_OUT, "panel_A_feature_stability.csv"))

# Permutation test on best classifier (1000 perms)
message(sprintf("  Permutation test on best classifier (%s)...", best_clf))

best_me <- switch(best_clf,
  Pre      = me_pre[common_subj, ],
  Post     = me_post[common_subj, ],
  DeltaME  = me_delta,
  Combined = me_avg
)
best_k_median <- as.integer(median(all_res[[best_clf]]$ks))
message(sprintf("  Using fixed k=%d (median) for permutations", best_k_median))

n_perm <- 1000
null_aucs <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  shuffled <- sample(true_labels)
  null_aucs[i] <- fast_loocv_auc(shuffled, best_me, best_k_median)
}

perm_pvalue <- (sum(null_aucs >= best_auc) + 1) / (n_perm + 1)
null_mean <- mean(null_aucs)
null_sd   <- sd(null_aucs)

message(sprintf("  %s AUC=%.3f, perm p=%.4f, null mean=%.3f +/-%.3f",
                best_clf, best_auc, perm_pvalue, null_mean, null_sd))

# Also run permutation on all 4 classifiers (100 perms each for speed)
message("  Quick permutation test on all classifiers (100 perms)...")
n_perm_quick <- 100
all_perm_p <- setNames(numeric(length(clf_names)), clf_names)

for (nm in clf_names) {
  me_nm <- switch(nm,
    Pre = me_pre[common_subj, ], Post = me_post[common_subj, ],
    DeltaME = me_delta, Combined = me_avg)
  k_nm <- as.integer(median(all_res[[nm]]$ks))
  obs_auc <- clf_summary |> filter(classifier == nm) |> pull(auc)
  quick_nulls <- numeric(n_perm_quick)
  for (j in seq_len(n_perm_quick)) {
    shuf <- sample(true_labels)
    quick_nulls[j] <- fast_loocv_auc(shuf, me_nm, k_nm)
  }
  all_perm_p[nm] <- (sum(quick_nulls >= obs_auc) + 1) / (n_perm_quick + 1)
}
all_perm_p[best_clf] <- perm_pvalue

clf_summary$perm_p <- all_perm_p[clf_summary$classifier]

write_csv(clf_summary, file.path(DAT_OUT, "panel_A_multi_classifier_auc.csv"))

perm_df <- tibble(
  best_classifier = best_clf,
  observed_auc    = best_auc,
  perm_pvalue     = perm_pvalue,
  null_auc_mean   = null_mean,
  null_auc_sd     = null_sd,
  null_auc_ci_lo  = null_mean - 1.96 * null_sd,
  null_auc_ci_hi  = null_mean + 1.96 * null_sd,
  n_permutations  = n_perm
)
write_csv(perm_df, file.path(DAT_OUT, "panel_A_permutation.csv"))


roc_all_df <- map_dfr(clf_names, function(nm) {
  probs <- all_res[[nm]]$probs
  roc_obj <- tryCatch(roc(true_labels, probs, quiet = TRUE),
                       error = function(e) NULL)
  if (is.null(roc_obj)) return(tibble())
  tibble(classifier = nm,
         fpr = 1 - roc_obj$specificities,
         tpr = roc_obj$sensitivities)
})
write_csv(roc_all_df, file.path(DAT_OUT, "panel_A_roc_curves.csv"))
