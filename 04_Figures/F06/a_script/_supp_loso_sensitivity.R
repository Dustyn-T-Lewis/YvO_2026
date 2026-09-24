#!/usr/bin/env Rscript
# F06 Supplementary: LOSO sensitivity for per-module univariate ROC AUCs
#
# Addresses eigengene-projection optimism (not module-definition circularity).
# For each cell the main figure draws, leave one subject out, refit the module's
# 1st PC on n-1 training subjects, project the held-out subject onto that PC,
# and compute an out-of-fold AUC.
#
# Sourced by F06_data.R after style.R and figure_supplement_helpers.R.

pacman::p_load(tidyverse, pROC)

BASE    <- "04_Figures/F06"
DAT_OUT <- file.path(BASE, "c_data", "loso_auc")
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

# Inputs
F05_SUPP  <- "04_Figures/F05/c_data/F05_data.xlsx"
stopifnot("F05 must run first: missing F05_data.xlsx" =
  file.exists(F05_SUPP))
datExpr   <- readRDS("04_Figures/F05/c_data/datExpr.rds")
mod_cols  <- readRDS("04_Figures/F05/c_data/module_colors.rds")
me_pre    <- readRDS("04_Figures/F05/c_data/me_pre.rds")
me_post   <- readRDS("04_Figures/F05/c_data/me_post.rds")
subj_age  <- read_sheet_df(F05_SUPP, "metadata_subj_age")
pheno     <- read_sheet_df(F05_SUPP, "metadata_pheno_wide")

# In-sample AUCs: prefer CSV (if S7_A_module_grid.R just ran), else xlsx
in_sample_csv  <- file.path(BASE, "c_data", "module_grid", "module_grid_summary.csv")
in_sample_xlsx <- file.path(BASE, "c_data", "F06_data.xlsx")
if (file.exists(in_sample_csv)) {
  in_sample <- read_csv(in_sample_csv, show_col_types = FALSE)
} else {
  in_sample <- read_sheet_df(in_sample_xlsx, "module_grid_summary")
}

sample_ids   <- rownames(datExpr)
subject_keys <- sub("_(Pre|Post)$", "", sample_ids)
common_subj  <- rownames(me_pre)
MODULES      <- c("turquoise","blue","brown","yellow","green","red","black","pink")

stopifnot(all(common_subj %in% subject_keys))

# The cells the main figure draws, not a separate ranking of the same 48.
# This table is read as the sensitivity check on the figure, so it has to
# cover the same cells.
source("04_Figures/F06/a_script/_panel_selection.R")
cells <- main_panel_cells(in_sample |> filter(!is.na(perm_p), !is.na(auc)))

message("Main-panel (row, module) pairs carried into LOSO:")
print(cells |> dplyr::select(row, module, auc, perm_p, q_bh))

# Outcome vectors (mirror _supp_module_grid.R)
age_bin <- ifelse(subj_age$age[match(common_subj, subj_age$subject_key)] == "Old", 1, 0)

pheno_s <- pheno[match(common_subj, pheno$subject_key), ]
ok_vl <- !is.na(pheno_s$delta_VL)
resid_vl <- rep(NA_real_, length(common_subj))
resid_vl[ok_vl] <- residuals(lm(pheno_s$delta_VL[ok_vl] ~ age_bin[ok_vl]))
vl_bin <- ifelse(resid_vl > median(resid_vl, na.rm = TRUE), 1L, 0L)

ok_lbm <- !is.na(pheno_s$LBM_Pre)
resid_lbm <- rep(NA_real_, length(common_subj))
resid_lbm[ok_lbm] <- residuals(lm(pheno_s$LBM_Pre[ok_lbm] ~ age_bin[ok_lbm]))
lbm_bin <- ifelse(resid_lbm > median(resid_lbm, na.rm = TRUE), 1L, 0L)

# LOSO projection helper
loso_me <- function(full_mat, train_rows, holdout_rows) {
  Xtr <- scale(full_mat[train_rows, , drop = FALSE])
  center <- attr(Xtr, "scaled:center")
  sc     <- attr(Xtr, "scaled:scale")
  sc[sc == 0 | !is.finite(sc)] <- 1
  sv <- svd(Xtr, nu = 0, nv = 1)
  v  <- sv$v[, 1]
  me_train <- as.numeric(Xtr %*% v)
  avg_train <- rowMeans(Xtr)
  if (cor(me_train, avg_train) < 0) v <- -v
  Xhold <- sweep(full_mat[holdout_rows, , drop = FALSE], 2, center, "-")
  Xhold <- sweep(Xhold, 2, sc, "/")
  as.numeric(Xhold %*% v)
}

stopifnot(length(mod_cols) == ncol(datExpr))
module_proteins <- function(module) colnames(datExpr)[mod_cols == module]

idx_pre   <- match(paste0(common_subj, "_Pre"),  rownames(datExpr))
idx_post  <- match(paste0(common_subj, "_Post"), rownames(datExpr))
stopifnot(!any(is.na(idx_pre)), !any(is.na(idx_post)))

X_pre_sub  <- datExpr[idx_pre,  , drop = FALSE]
X_post_sub <- datExpr[idx_post, , drop = FALSE]
X_comb_sub <- (X_pre_sub + X_post_sub) / 2
rownames(X_pre_sub) <- rownames(X_post_sub) <- rownames(X_comb_sub) <- common_subj

# Per-pair LOSO
results <- list()

set.seed(42)
for (i in seq_len(nrow(cells))) {
  row <- cells$row[i]; mod <- cells$module[i]
  mp  <- module_proteins(mod)
  if (length(mp) < 2) {
    message(sprintf("[skip] %s / %s -- only %d proteins in module", row, mod, length(mp)))
    next
  }

  if (row == "Age") {
    X_sub <- X_comb_sub[, mp, drop = FALSE]; y <- age_bin; ok <- rep(TRUE, length(y))
    paired <- FALSE
  } else if (row == "\u0394VL-responder") {
    X_sub <- X_pre_sub[, mp, drop = FALSE]; y <- vl_bin; ok <- !is.na(y)
    paired <- FALSE
  } else if (row == "LBM (High vs Low)") {
    X_sub <- X_pre_sub[, mp, drop = FALSE]; y <- lbm_bin; ok <- !is.na(y)
    paired <- FALSE
  } else if (row == "Pre vs Post") {
    paired <- TRUE; age_mask <- rep(TRUE, length(common_subj))
  } else if (row == "Pre vs Post (Young)") {
    paired <- TRUE; age_mask <- age_bin == 0
  } else if (row == "Pre vs Post (Old)") {
    paired <- TRUE; age_mask <- age_bin == 1
  } else {
    message(sprintf("[skip] %s -- unrecognised row", row)); next
  }

  if (paired) {
    subj_use <- common_subj[age_mask]
    Xp_full  <- X_pre_sub[subj_use,  mp, drop = FALSE]
    Xq_full  <- X_post_sub[subj_use, mp, drop = FALSE]
    X_pool <- rbind(Xp_full, Xq_full)
    tp_lab <- c(rep(0L, nrow(Xp_full)), rep(1L, nrow(Xq_full)))
    sid    <- c(rownames(Xp_full), rownames(Xq_full))

    preds <- rep(NA_real_, length(tp_lab))
    for (s in unique(sid)) {
      hold <- which(sid == s); train <- which(sid != s)
      preds[hold] <- loso_me(X_pool, train, hold)
    }
    preds_ok <- !is.na(preds)
    r <- suppressMessages(roc(tp_lab[preds_ok], preds[preds_ok],
                               quiet = TRUE, direction = "auto"))
    ci <- as.numeric(ci.auc(r))
    results[[length(results)+1]] <- tibble(
      row = row, module = mod,
      n_subj  = length(unique(sid)),
      n_obs   = sum(preds_ok),
      auc_insample = cells$auc[i],
      auc_loso     = as.numeric(auc(r)),
      ci_lo_loso   = ci[1], ci_hi_loso = ci[3],
      drop = cells$auc[i] - as.numeric(auc(r)))
  } else {
    idx <- which(ok)
    preds <- rep(NA_real_, length(idx))
    y_vec <- y[idx]
    Xu    <- X_sub[idx, , drop = FALSE]
    for (k in seq_along(idx)) {
      hold <- k; train <- setdiff(seq_along(idx), k)
      preds[k] <- loso_me(Xu, train, hold)
    }
    preds_ok <- !is.na(preds)
    if (length(unique(y_vec[preds_ok])) < 2) {
      message(sprintf("[skip] %s / %s -- held-out labels collapse", row, mod)); next
    }
    r <- suppressMessages(roc(y_vec[preds_ok], preds[preds_ok],
                               quiet = TRUE, direction = "auto"))
    ci <- as.numeric(ci.auc(r))
    results[[length(results)+1]] <- tibble(
      row = row, module = mod,
      n_subj  = sum(preds_ok),
      n_obs   = sum(preds_ok),
      auc_insample = cells$auc[i],
      auc_loso     = as.numeric(auc(r)),
      ci_lo_loso   = ci[1], ci_hi_loso = ci[3],
      drop = cells$auc[i] - as.numeric(auc(r)))
  }
}

loso_df <- bind_rows(results) |>
  mutate(row_module = paste(module, row, sep = " | ")) |>
  arrange(desc(auc_insample))

write_csv(loso_df, file.path(DAT_OUT, "loso_auc_summary.csv"))

message("LOSO sensitivity: in-sample vs LOSO AUC")
print(loso_df |> dplyr::select(row, module, n_subj, n_obs,
                                auc_insample, auc_loso, drop))
message(sprintf("Median in-sample AUC: %.3f", median(loso_df$auc_insample)))
message(sprintf("Median LOSO AUC:      %.3f", median(loso_df$auc_loso)))
message(sprintf("Median drop:          %.3f", median(loso_df$drop)))

message("Wrote: ", file.path(DAT_OUT, "loso_auc_summary.csv"))
