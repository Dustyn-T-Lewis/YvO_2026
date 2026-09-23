#!/usr/bin/env Rscript
# Build the ranked-list fGSEA feature cache from the DEP t-statistics.
# Sourced by any figure that needs it; rebuilt automatically when the DEP
# results are newer than the cache. It used to be skipped whenever the file
# existed, which let the cache silently predate the DEP fit.

setwd(here::here())

CACHE_PATH <- "04_Figures/shared/fgsea_tstat_all_v2.csv"
STAGE3_CSV <- "03_DEP/c_data/03_combined_results.csv"

# Gene-set size bounds as reported in Methods. One definition, used for both the
# collection and the fgsea call, so the two cannot drift apart.
MIN_SIZE <- 15L
MAX_SIZE <- 500L
PADJ_CUTOFF <- 0.05

stopifnot(
  "Stage 03 combined_results.csv missing — run 03_DEP/a_script/01_run_dep.R" =
    file.exists(STAGE3_CSV)
)

fresh <- file.exists(CACHE_PATH) &&
  file.mtime(CACHE_PATH) >= file.mtime(STAGE3_CSV)
if (fresh) {
  message(sprintf(
    "fGSEA cache present (%s) — up to date with the DEP results",
    basename(CACHE_PATH)
  ))
} else {
  message("Rebuilding fGSEA cache from current Stage 03 t-statistics...")

  source("04_Figures/shared/pathway_utils.R")
  source("04_Figures/shared/enrichment_dedup.R")
  pacman::p_load(dplyr, readr, tibble)

  CTRS <- c("Aging", "Training_Young", "Training_Old", "Interaction")
  dep <- read_csv(STAGE3_CSV, show_col_types = FALSE)

  pw_list <- build_pathway_collection(
    min_size = MIN_SIZE, max_size = MAX_SIZE,
    include_goslim = TRUE,
    exclude_variants = FALSE
  )

  # Methods states Benjamini-Hochberg correction per database, so each
  # collection is tested as its own family rather than pooled.
  pw_by_db <- split(pw_list, classify_database(names(pw_list)))

  set.seed(42)
  all_results <- list()
  for (ctr in CTRS) {
    message(sprintf("--- %s ---", ctr))
    ranks <- setNames(dep[[paste0("t_", ctr)]], dep$gene)
    ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

    per_db <- lapply(pw_by_db, \(pw) run_fgsea_deduplicated(
      ranks     = ranks,
      pathways  = pw,
      em_cutoff = 1, # no-op here -- the cross-database em_dedup_report() call below does the real collapse
      nperm     = 10000,
      min_size  = MIN_SIZE,
      max_size  = MAX_SIZE
    ))
    # em_cutoff = 1 above makes run_fgsea_deduplicated()'s own flagging a
    # no-op, but it still stamps the three dedup columns onto every row --
    # drop them so the real, cross-database flags below don't collide with
    # them in the join.
    res <- bind_rows(per_db) |>
      select(-any_of(c("dedup_status", "merged_into", "overlap_jaccard")))

    # Redundancy is flagged, never dropped: pathway-count panels report every
    # significant set, rings show one representative per redundant cluster.
    sig <- as.data.frame(res[!is.na(res$padj) & res$padj < PADJ_CUTOFF, ])
    flags <- em_dedup_report(sig, pw_list)[, c(
      "pathway", "dedup_status", "merged_into", "overlap_jaccard"
    )]

    all_results[[ctr]] <- res |>
      left_join(flags, by = "pathway") |>
      mutate(contrast = ctr)
  }

  cache_df <- bind_rows(all_results) |>
    mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))
  cache_df <- cache_df[, c(
    "pathway", "pval", "padj", "log2err", "ES", "NES",
    "size", "leadingEdge", "database", "contrast",
    "dedup_status", "merged_into", "overlap_jaccard"
  )]

  write_csv(cache_df, CACHE_PATH)
  message(sprintf(
    "Wrote %d rows to %s (%d significant, %d kept after EnrichmentMap collapse)",
    nrow(cache_df), CACHE_PATH,
    sum(cache_df$padj < PADJ_CUTOFF, na.rm = TRUE),
    sum(cache_df$dedup_status == "kept", na.rm = TRUE)
  ))
}
