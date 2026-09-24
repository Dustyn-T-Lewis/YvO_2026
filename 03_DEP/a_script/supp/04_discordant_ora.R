# 04_discordant_ora: R2.2's second half: do the whole-proteome discordant
# proteins group into biological pathways?
#
# The reviewer's ~40% is the proteins whose training logFCs disagree in sign
# between age groups, counted among those with an estimate in both.
# _supp_concordance_magnitude.R showed sign disagreement is mostly a
# small-effect phenomenon: three quarters involve a protein whose smaller
# effect sits below |log2FC| = 0.1. Pre-specified 2026-08-18: the primary test
# is ORA on the discordant proteins with |log2FC| >= 0.1 in BOTH groups (a
# real effect in opposite directions); the full discordant set and the
# direction-split subsets are reported alongside. Universe is every protein
# with an estimate in both groups; enrichment machinery
# and dedup parameters match the figure pipeline (fora + EnrichmentMap-style
# collapse, Jaccard 0.5, set size 10-500).

pacman::p_load(withr, readr, dplyr, tidyr, tibble, purrr)

withr::local_dir(here::here())
source("04_Figures/shared/pathway_utils.R")

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

MIN_EFFECT <- 0.1

both <- dep |>
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) |>
  mutate(
    discordant = sign(logFC_Training_Young) != sign(logFC_Training_Old),
    real_both = pmin(abs(logFC_Training_Young), abs(logFC_Training_Old)) >= MIN_EFFECT,
    pattern = if_else(logFC_Training_Young > 0, "Y up / O down", "Y down / O up")
  )

universe <- unique(both$gene)
disc <- filter(both, discordant)

message(sprintf(
  "%d of %d discordant (%.1f%%); %d with |log2FC| >= %.1f in both groups",
  nrow(disc), nrow(both), 100 * nrow(disc) / nrow(both),
  sum(disc$real_both), MIN_EFFECT
))

pw_list <- build_pathway_collection(
  min_size = 10, max_size = 500,
  include_goslim = FALSE, exclude_variants = TRUE
)

# run_ora_deduplicated errors when nothing reaches significance; an empty
# result is a finding here, not a failure.
ora_on <- function(genes, set_label) {
  tryCatch(
    run_ora_deduplicated(
      genes = genes, universe = universe, pathways = pw_list,
      em_cutoff = 0.5, min_size = 10, max_size = 500, padj_cutoff = 1
    ) |>
      mutate(set = set_label, n_genes = length(genes)),
    error = function(e) {
      message(set_label, ": no enriched sets")
      tibble(set = set_label, n_genes = length(genes))
    }
  )
}

set.seed(42)
disc_ora <- bind_rows(
  ora_on(disc$gene[disc$real_both], sprintf("discordant, both |logFC| >= %.1f", MIN_EFFECT)),
  ora_on(disc$gene, sprintf("discordant, all %d", nrow(disc))),
  ora_on(
    disc$gene[disc$real_both & disc$pattern == "Y up / O down"],
    "real discordant, Y up / O down"
  ),
  ora_on(
    disc$gene[disc$real_both & disc$pattern == "Y down / O up"],
    "real discordant, Y down / O up"
  )
)

write_csv(disc_ora, "03_DEP/c_data/04_discordant_ora.csv")

if ("dedup_status" %in% names(disc_ora)) {
  disc_ora |>
    filter(dedup_status == "kept", padj < 0.05) |>
    summarise(n_terms = n(), top = pathway[which.min(padj)], .by = set) |>
    as.data.frame() |>
    print()
} else {
  message("No set returned any enriched pathway.")
}
