# Redundancy collapsing for pathway enrichment results.
#
# Greedy in padj order, so the most significant member of a redundant cluster
# becomes its representative. Sets are flagged, never dropped, so one
# computation serves both views: stacked bars count every significant set,
# rings show representatives only. Similarity uses the full gene sets, so
# redundancy reflects database structure and is contrast-independent.
#
# Jaccard at 0.5, within database, as in the submitted paper. EnrichmentMap's
# own coefficient (Merico 2010 PMID 21085593; Reimand 2019 PMID 30664679 Box 1)
# averages Jaccard with the Szymkiewicz-Simpson overlap coefficient and
# thresholds at 0.375, which collapses roughly twice as many terms.

jaccard <- function(genes_a, genes_b) {
  inter <- length(intersect(genes_a, genes_b))
  if (inter == 0L) {
    return(0)
  }
  inter / (length(genes_a) + length(genes_b) - inter)
}

# Names of the pathways that survive collapsing, most significant first.
collapse_redundant <- function(results, pathways, cutoff) {
  results <- results[order(results$padj), ]
  kept_genes <- list()
  keep <- logical(nrow(results))
  for (i in seq_len(nrow(results))) {
    genes <- pathways[[results$pathway[i]]]
    if (is.null(genes)) {
      keep[i] <- TRUE
      next
    }
    dup <- vapply(kept_genes, \(k) jaccard(genes, k) >= cutoff, logical(1))
    if (!any(dup)) {
      keep[i] <- TRUE
      kept_genes[[length(kept_genes) + 1L]] <- genes
    }
  }
  results$pathway[keep]
}

# Flags each set as representative or redundant; for redundant ones, names the
# representative it overlaps most and the Jaccard to it.
em_dedup_report <- function(results, pathways, cutoff = 0.5) {
  if (nrow(results) == 0) {
    results$dedup_status <- character(0)
    results$merged_into <- character(0)
    results$overlap_jaccard <- numeric(0)
    return(results)
  }
  by_db <- if ("database" %in% names(results)) results$database else "all"
  kept <- unlist(
    lapply(split(results, by_db), collapse_redundant,
      pathways = pathways, cutoff = cutoff
    ),
    use.names = FALSE
  )
  redundant <- setdiff(results$pathway, kept)

  results$dedup_status <- ifelse(
    results$pathway %in% redundant, "redundant", "kept"
  )
  results$merged_into <- NA_character_
  results$overlap_jaccard <- NA_real_
  for (p in redundant) {
    genes <- pathways[[p]]
    if (is.null(genes)) next
    sims <- vapply(kept, \(b) jaccard(genes, pathways[[b]]), numeric(1))
    best <- which.max(sims)
    results$merged_into[results$pathway == p] <- kept[best]
    results$overlap_jaccard[results$pathway == p] <- sims[best]
  }
  results[order(results$padj), , drop = FALSE]
}
