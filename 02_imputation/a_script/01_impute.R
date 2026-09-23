#!/usr/bin/env Rscript
# Stage 02: Impute
# MAR/MNAR classification (k-means on intensity x missingness) → missForest
#
# Outputs:
#   c_data/02_imputation.xlsx     5-sheet supplement (+benchmark if available)
#   c_data/01_DAList_imputed.rds  proteoDA object with imputation annotations
#   c_data/00_report_intermediates.rds diagnostic data for reports + F00

withr::local_dir(here::here())

pacman::p_load(missForest, dplyr, tibble, openxlsx)
source("04_Figures/shared/devices.R")
source("04_Figures/shared/supplement_overview.R")

set.seed(42)

DAT <- "02_imputation/c_data"
dir.create(DAT, showWarnings = FALSE, recursive = TRUE)

MISS_UNRELIABLE <- 50 # % missing above which imputation is flagged unreliable

# Palettes (passed to reports via intermediates)

PAL_GT <- c(
  Young_Pre = scales::alpha("#4393C3", 0.5), Young_Post = "#4393C3",
  Old_Pre = scales::alpha("#D6604D", 0.5), Old_Post = "#D6604D"
)
PAL_MAR <- c(MAR = "#4393C3", MNAR = "#D6604D")
PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#4393C3", MNAR = "#D6604D")

# Load from Stage 01
# Read numeric matrix from CSV (text-serialized) for cross-pipeline
# reproducibility; RDS binary doubles differ at machine-epsilon from CSV
# round-trip, which shifts missForest tree splits.

df <- readr::read_csv("01_normalization/c_data/02_normalized.csv",
  show_col_types = FALSE
)
ann <- df |> select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- df$gene

if (any(duplicated(df$gene))) {
  warning("Duplicate gene names; using uniprot_id as rownames")
  rownames(mat) <- df$uniprot_id
}

dal_norm <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
meta <- as_tibble(dal_norm$metadata) |>
  select(Col_ID, Group, Timepoint, Group_Time)

stopifnot(setequal(meta$Col_ID, colnames(mat)))
message(sprintf("Loaded: %d proteins x %d samples", nrow(mat), ncol(mat)))

# Missingness profiling

prot_miss <- rowSums(is.na(mat))
prot_pct <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss <- round(sum(is.na(mat)) / length(mat) * 100, 2)

miss_by_group <- sapply(unique(meta$Group_Time), \(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

# MAR/MNAR classification
#
# One rule, k-means on (mean intensity, % missing). MNAR is the cluster with
# the lower mean intensity: a protein missing because it sits near the
# detection limit. This replaced a "3-method consensus" on 2026-09-08 in which
# two of the three votes were forced median splits and so labelled exactly
# 49.96% of incomplete proteins MNAR whatever the data showed. The logistic
# classifier fitted P(missing | intensity) with one continuous predictor and
# thresholded at the median of its own fitted values, which is algebraically a
# median split on intensity -- identical on all 1151 incomplete proteins. The
# left-tail classifier thresholded at its own median too. Only this one
# responds to structure in the data.

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
inc_mean <- obs_means[has_na]
inc_pct <- prot_pct[has_na]

set.seed(42)
km <- kmeans(scale(cbind(inc_mean, inc_pct)), centers = 2, nstart = 25)
km_mnar <- km$cluster == which.min(tapply(inc_mean, km$cluster, mean))

miss_class <- tibble(
  gene           = rownames(mat),
  n_miss         = prot_miss,
  pct_miss       = prot_pct,
  mean_intensity = obs_means,
  mnar_cluster   = NA_integer_
)
miss_class$mnar_cluster[has_na] <- as.integer(km_mnar)

miss_class <- miss_class |>
  mutate(
    classification = case_when(
      n_miss == 0 ~ "Complete",
      mnar_cluster == 1L ~ "MNAR",
      TRUE ~ "MAR"
    ),
    imputation_reliable = classification == "Complete" | pct_miss < MISS_UNRELIABLE
  )

mnar_genes <- miss_class$gene[miss_class$classification == "MNAR"]

# Group-stratified Fisher test per MNAR protein
set.seed(42)
group_miss_pval <- setNames(rep(NA_real_, nrow(miss_class)), miss_class$gene)
for (g in mnar_genes) {
  ct <- sapply(unique(meta$Group_Time), \(gt) {
    cols <- meta$Col_ID[meta$Group_Time == gt]
    c(missing = sum(is.na(mat[g, cols])), observed = sum(!is.na(mat[g, cols])))
  })
  group_miss_pval[g] <- tryCatch(
    fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
    error = \(e) NA_real_
  )
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$gene]

n_mar <- sum(miss_class$classification == "MAR")
n_mnar <- length(mnar_genes)
n_comp <- sum(miss_class$classification == "Complete")
mar_vals <- sum(miss_class$n_miss[miss_class$classification == "MAR"])
mnar_vals <- sum(miss_class$n_miss[miss_class$classification == "MNAR"])
total_vals <- mar_vals + mnar_vals

message(sprintf("Classification: MAR %d | MNAR %d | Complete %d", n_mar, n_mnar, n_comp))

# missForest imputation

message("Imputing with missForest...")
gene_order <- order(rownames(mat))
mat <- mat[gene_order, ]
ann <- ann[gene_order, ]
set.seed(42)
mf <- missForest::missForest(t(mat), maxiter = 10, ntree = 100, verbose = TRUE)
mat_imp <- t(mf$ximp)
rownames(mat_imp) <- rownames(mat)
colnames(mat_imp) <- colnames(mat)
stopifnot(sum(is.na(mat_imp)) == 0)
oob <- as.numeric(mf$OOBerror[1])
message(sprintf("OOB error: %.4f", oob))

# MNAR audit

was_na <- is.na(mat)

mnar_audit <- tibble(
  gene = mnar_genes,
  pre_mean = rowMeans(mat[mnar_genes, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[mnar_genes, ]),
  pre_sd = apply(mat[mnar_genes, ], 1, sd, na.rm = TRUE),
  pct_miss = prot_pct[mnar_genes],
  shift = post_mean - pre_mean,
  effect_d = (post_mean - pre_mean) / pre_sd,
  imputation_reliable = prot_pct[mnar_genes] < MISS_UNRELIABLE
)

# Build xlsx


imp_df <- bind_cols(ann, as_tibble(mat_imp))
mask_df <- bind_cols(tibble(gene = rownames(was_na)), as_tibble(was_na + 0L))
summary_df <- tibble(
  metric = c(
    "n_proteins", "n_samples", "pct_missing", "n_complete",
    "n_mar_proteins", "n_mnar_proteins", "n_mar_values", "n_mnar_values",
    "method", "oob_error"
  ),
  value = c(
    nrow(mat), ncol(mat), pct_miss, n_comp, n_mar, n_mnar,
    mar_vals, mnar_vals, "missForest", round(oob, 4)
  )
)

wb <- createWorkbook()
write_sheet(wb, "imputed_matrix", imp_df)
write_sheet(wb, "mar_mnar_classification", as.data.frame(miss_class))
write_sheet(wb, "imputation_mask", mask_df)
write_sheet(wb, "mnar_audit", as.data.frame(mnar_audit))
write_sheet(wb, "imputation_summary", summary_df)

bm_path <- file.path(DAT, "benchmark", "04_composite_ranking.csv")
if (file.exists(bm_path)) {
  bm <- readr::read_csv(bm_path, show_col_types = FALSE)
  write_sheet(wb, "benchmark_ranking", as.data.frame(bm))
}

add_overview(
  wb,
  title = "S9 Table \u2014 imputation",
  description = paste(
    "The imputed matrix the models were fitted on, which values were imputed,",
    "and the benchmark that chose the method."
  ),
  entries = c(
    imputed_matrix          = "The matrix after missForest imputation, one row per protein and one column per sample.",
    mar_mnar_classification = "Whether each protein's missing values look random or abundance-dependent.",
    imputation_mask         = "A 1 for every value that was imputed and a 0 for every value that was measured.",
    mnar_audit              = "Checks on the abundance-dependent proteins, which imputation can bias.",
    imputation_summary      = "How many values were imputed, and the method's out-of-bag error.",
    benchmark_ranking       = "Sixteen imputation methods scored on reconstruction, downstream effect and stability, against an unimputed reference row."
  )
)
saveWorkbook(wb, file.path(DAT, "02_imputation.xlsx"), overwrite = TRUE)

# CSVs for benchmark infrastructure and downstream figures
readr::write_csv(imp_df, file.path(DAT, "01_imputed.csv"))
readr::write_csv(as.data.frame(miss_class), file.path(DAT, "02_mar_mnar_classification.csv"))

# Save R objects

mat_imp_uid <- mat_imp
rownames(mat_imp_uid) <- ann$uniprot_id
dal <- dal_norm
dal$data <- mat_imp_uid
n_ann <- nrow(dal$annotation)
dal$annotation <- merge(
  dal$annotation,
  miss_class |> select(gene, n_miss, pct_miss,
    miss_classification = classification,
    imputation_reliable
  ),
  by = "gene", all.x = TRUE, sort = FALSE
)
stopifnot(nrow(dal$annotation) == n_ann)
# Re-align $annotation rows to $data row order. mat_imp was reordered by
# gene_order for missForest determinism; merge() preserves left-frame order.
# Without this match() step the saved DAList has the same set of proteins
# in $data and $annotation but at different row positions.
dal$annotation <- dal$annotation[
  match(rownames(dal$data), dal$annotation$uniprot_id), ,
  drop = FALSE
]
rownames(dal$annotation) <- dal$annotation$uniprot_id
stopifnot(identical(rownames(dal$data), dal$annotation$uniprot_id))
saveRDS(dal, file.path(DAT, "01_DAList_imputed.rds"))

saveRDS(list(
  mat = mat, mat_imp = mat_imp, was_na = was_na, ann = ann, meta = meta,
  miss_class = miss_class, miss_by_group = miss_by_group,
  prot_pct = prot_pct, pct_miss = pct_miss,
  mnar_genes = mnar_genes, mnar_audit = mnar_audit,
  n_mar_prots = n_mar, n_mnar_prots = n_mnar,
  mar_miss_vals = mar_vals, mnar_miss_vals = mnar_vals,
  total_miss_vals = total_vals, oob_error = oob,
  classification_method = "k-means on mean intensity x percent missing",
  PAL_GT = PAL_GT, PAL_MAR = PAL_MAR, PAL_CLASS = PAL_CLASS
), file.path(DAT, "00_report_intermediates.rds"))

message(sprintf(
  "Done: %d proteins x %d samples | OOB=%.4f",
  nrow(mat_imp), ncol(mat_imp), oob
))
