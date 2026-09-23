# Unified pathway enrichment utilities
# MSigDB Hallmark (H), Canonical Pathways (C2:CP), GO:BP (C5:GO:BP)
#
# Exports:
#   build_pathway_collection()      assemble GO:BP + Reactome + Hallmark + KEGG
#   run_fgsea_deduplicated()        fgsea + collapsePathways for one contrast
#   run_enrichment_pipeline()       fGSEA across all contrasts + databases
#   run_ora_deduplicated()          over-representation, redundancy flagged
#   classify_database() / classify_pathway_func()  category labels for plotting

source(here::here("04_Figures", "shared", "enrichment_dedup.R"))

deduplicate_enrichment_flat <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) {
    return(results)
  }

  results <- results[order(results$padj), ]
  kept_sets <- list()
  keep_mask <- logical(nrow(results))

  for (i in seq_len(nrow(results))) {
    pw_name <- results$pathway[i]
    pw_genes <- pathways[[pw_name]]
    if (is.null(pw_genes)) {
      keep_mask[i] <- TRUE
      next
    }

    is_redundant <- FALSE
    for (j in seq_along(kept_sets)) {
      inter <- length(intersect(pw_genes, kept_sets[[j]]))
      union <- length(union(pw_genes, kept_sets[[j]]))
      if (union > 0 && (inter / union) > jaccard_cutoff) {
        is_redundant <- TRUE
        break
      }
    }

    if (!is_redundant) {
      keep_mask[i] <- TRUE
      kept_sets[[length(kept_sets) + 1]] <- pw_genes
    }
  }

  results[keep_mask, ]
}

# Database-stratified dedup: within-db first, then merge survivors by padj.
# Falls back to flat dedup if no 'database' column.
deduplicate_enrichment <- function(results, pathways, jaccard_cutoff = 0.5) {
  if (nrow(results) == 0) {
    return(results)
  }

  if (!"database" %in% names(results)) {
    return(deduplicate_enrichment_flat(results, pathways, jaccard_cutoff))
  }

  dbs <- unique(results$database)
  within_dedup <- list()
  for (db in dbs) {
    db_rows <- results[results$database == db, ]
    within_dedup[[db]] <- deduplicate_enrichment_flat(
      db_rows, pathways,
      jaccard_cutoff
    )
  }

  survivors <- do.call(rbind, within_dedup)
  survivors[order(survivors$padj), ]
}


build_pathway_collection <- function(species = "Homo sapiens",
                                     min_size = 10, max_size = 500,
                                     include_goslim = TRUE,
                                     exclude_variants = FALSE) {
  requireNamespace("msigdbr", quietly = TRUE)
  stopifnot(
    "msigdbr release differs from the pinned 26.1.0; every enrichment result will move" =
      as.character(packageVersion("msigdbr")) == "26.1.0"
  )

  hallmark <- msigdbr::msigdbr(species = species, collection = "H")
  kegg <- msigdbr::msigdbr(
    species = species, collection = "C2",
    subcollection = "CP:KEGG_MEDICUS"
  )
  reactome <- msigdbr::msigdbr(
    species = species, collection = "C2",
    subcollection = "CP:REACTOME"
  )
  gobp <- msigdbr::msigdbr(
    species = species, collection = "C5",
    subcollection = "GO:BP"
  )

  # Applied to every collection, so the exclusion Methods reports holds without
  # qualification. Pattern matches Mito.
  disease_pat <- paste0(
    "DISEASE|CANCER|TUMOR|CARCINOMA|LEUKEMIA|LYMPHOMA|MELANOMA|GLIOMA|",
    "HEPATITIS|HIV|INFECTION|INFECTIOUS|VIRAL|VIRUS|INFLUENZA|SARS|HCMV|",
    "MEASLES|DENGUE|BACTERIAL|LISTERIA|LEISHMANIA|PARASIT"
  )
  drop_disease <- function(df) {
    df[!grepl(disease_pat, df$gs_name, ignore.case = TRUE), ]
  }
  hallmark <- drop_disease(hallmark)
  kegg <- drop_disease(kegg)
  reactome <- drop_disease(reactome)
  gobp <- drop_disease(gobp)

  if (exclude_variants) {
    kegg <- kegg[!grepl("_VARIANT_", kegg$gs_name), ]
  }

  cols <- c("gs_name", "gene_symbol")
  sets_list <- list(hallmark[, cols], kegg[, cols], reactome[, cols], gobp[, cols])
  dbs <- c("H", "KEGG", "Reactome", "GO:BP")
  all_sets <- do.call(rbind, sets_list)

  pw_list <- split(all_sets$gene_symbol, all_sets$gs_name)
  pw_list <- lapply(pw_list, unique)

  if (include_goslim) {
    goslim_sets <- build_goslim_gene_sets(
      species = species, min_size = min_size, max_size = max_size
    )
    goslim_sets <- goslim_sets[
      !grepl(disease_pat, names(goslim_sets), ignore.case = TRUE)
    ]
    pw_list <- c(pw_list, goslim_sets)
    dbs <- c(dbs, "GO Slim")
  }

  sizes <- vapply(pw_list, length, integer(1))
  pw_list <- pw_list[sizes >= min_size & sizes <= max_size]

  message(sprintf(
    "Pathway collection: %d sets (%s), size %d-%d",
    length(pw_list), paste(dbs, collapse = " + "),
    min_size, max_size
  ))
  pw_list
}


# Official GO Consortium generic slim, checked in from
# https://current.geneontology.org/ontology/subsets/goslim_generic.obo
# (data-version go/releases/2026-07-26/subsets/goslim_generic.owl, verified
# against the live file 2026-08-20). Every non-obsolete biological_process
# term in that file, nothing hand-picked or trimmed.
parse_goslim_bp_terms <- function(obo_path) {
  stopifnot("GO Slim OBO file missing" = file.exists(obo_path))
  lines <- readLines(obo_path)

  # Stanza boundaries are any bracketed header ([Term], [Typedef], ...), not
  # just [Term] -- otherwise the last [Term] block swallows every [Typedef]
  # stanza appended after it.
  stanza_start <- grep("^\\[.*\\]$", lines)
  stanza_end <- c(stanza_start[-1] - 1, length(lines))
  term_idx <- which(lines[stanza_start] == "[Term]")

  ids <- vapply(term_idx, function(i) {
    block <- lines[stanza_start[i]:stanza_end[i]]
    ns <- sub("^namespace: ", "", block[startsWith(block, "namespace: ")])
    obsolete <- any(startsWith(block, "is_obsolete: true"))
    if (length(ns) && ns == "biological_process" && !obsolete) {
      sub("^id: ", "", block[startsWith(block, "id: ")][1])
    } else {
      NA_character_
    }
  }, character(1))

  unique(ids[!is.na(ids)])
}

build_goslim_gene_sets <- function(species = "Homo sapiens",
                                   min_size = 10, max_size = 500,
                                   obo_path = here::here("04_Figures", "shared", "goslim_generic.obo")) {
  requireNamespace("GO.db", quietly = TRUE)
  requireNamespace("org.Hs.eg.db", quietly = TRUE)
  requireNamespace("AnnotationDbi", quietly = TRUE)

  bp_slim <- parse_goslim_bp_terms(obo_path)

  # Get all descendant GO terms for each slim term
  offspring <- as.list(GO.db::GOBPOFFSPRING)

  # Map all BP GO terms -> gene symbols via org.Hs.eg.db
  suppressMessages({
    go_genes <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "GO"),
      keytype = "GO",
      columns = c("SYMBOL", "ONTOLOGY")
    )
  })
  go_bp_genes <- go_genes[!is.na(go_genes$ONTOLOGY) & go_genes$ONTOLOGY == "BP", ]
  go_to_symbols <- split(go_bp_genes$SYMBOL, go_bp_genes$GO)

  # Build gene sets: each slim term + all its descendants
  goslim_sets <- list()
  slim_names <- vapply(bp_slim, function(id) {
    tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[id]]),
      error = function(e) NA_character_
    )
  }, character(1))

  for (i in seq_along(bp_slim)) {
    go_id <- bp_slim[i]
    go_term <- slim_names[i]
    if (is.na(go_term)) next

    # Collect genes from this term + all offspring
    all_terms <- go_id
    desc <- offspring[[go_id]]
    if (!is.null(desc)) all_terms <- c(all_terms, desc)

    genes <- unique(unlist(go_to_symbols[intersect(all_terms, names(go_to_symbols))],
      use.names = FALSE
    ))
    genes <- genes[!is.na(genes)]

    if (length(genes) >= min_size && length(genes) <= max_size) {
      set_name <- paste0("GOSLIM_", toupper(gsub(" ", "_", go_term)))
      goslim_sets[[set_name]] <- genes
    }
  }

  message(sprintf(
    "GO Slim: %d/%d terms passed size filter (%d-%d)",
    length(goslim_sets), length(bp_slim), min_size, max_size
  ))
  goslim_sets
}



# Redundancy is flagged, not dropped -- matches run_ora_deduplicated(), so a
# count like "19 significant" and a display of "representatives only" come
# from the same computation instead of two different redundancy rules. Uses
# Jaccard at 0.5, within database, matching the submitted paper.
run_fgsea_deduplicated <- function(ranks, pathways, em_cutoff = 0.5,
                                   nperm = 10000, min_size = 15,
                                   max_size = 500) {
  requireNamespace("fgsea", quietly = TRUE)

  # Protein groups can share a gene symbol; fgsea needs one stat per gene.
  # Keep the strongest signal.
  if (anyDuplicated(names(ranks))) {
    ranks <- ranks[order(-abs(ranks))]
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)
  }

  set.seed(42) # fgseaMultilevel is permutation-based; seed for reproducible NES/p
  res <- fgsea::fgseaMultilevel(
    pathways    = pathways,
    stats       = ranks,
    minSize     = min_size,
    maxSize     = max_size,
    nPermSimple = nperm,
    eps         = 0
  )
  res <- as.data.frame(res)

  res$database <- classify_database(res$pathway)

  res <- tibble::as_tibble(res)

  keep_cols <- c(
    "pathway", "padj", "NES", "size", "leadingEdge",
    "database", "pval", "ES", "log2err"
  )
  res <- res[, intersect(keep_cols, names(res))]

  sig <- res[!is.na(res$padj) & res$padj < 0.05, ]
  nonsig <- res[is.na(res$padj) | res$padj >= 0.05, ]

  sig_flagged <- em_dedup_report(sig, pathways, em_cutoff)
  nonsig[c("dedup_status", "merged_into", "overlap_jaccard")] <-
    list(NA_character_, NA_character_, NA_real_)

  n_kept <- sum(sig_flagged$dedup_status == "kept", na.rm = TRUE)
  message(sprintf(
    "fGSEA dedup: %d sig -> %d kept (removed %d, %.1f%%)",
    nrow(sig), n_kept, nrow(sig) - n_kept,
    if (nrow(sig) > 0) round(100 * (nrow(sig) - n_kept) / nrow(sig), 1) else 0
  ))

  rbind(sig_flagged, nonsig)
}


run_enrichment_pipeline <- function(stats_list, pw_list,
                                    jaccard_cutoff = 0.35,
                                    nperm = 10000,
                                    min_size = 15, max_size = 500,
                                    padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  all_results <- list()

  for (ctr in names(stats_list)) {
    message(sprintf("\n--- %s ---", ctr))
    ranks <- stats_list[[ctr]]
    if (anyDuplicated(names(ranks))) {
      # protein groups can share a gene symbol; keep the strongest stat
      ranks <- ranks[order(-abs(ranks))]
      ranks <- ranks[!duplicated(names(ranks))]
      ranks <- sort(ranks, decreasing = TRUE)
    }

    set.seed(42) # fgseaMultilevel is permutation-based; seed for reproducible NES/p
    res_dt <- fgsea::fgseaMultilevel(
      pathways    = pw_list,
      stats       = ranks,
      minSize     = min_size,
      maxSize     = max_size,
      nPermSimple = nperm,
      eps         = 0
    )

    sig_dt <- res_dt[!is.na(res_dt$padj) & res_dt$padj < padj_cutoff, ]
    if (nrow(sig_dt) > 0) {
      collapsed <- fgsea::collapsePathways(
        fgseaRes     = sig_dt,
        pathways     = pw_list,
        stats        = ranks
      )
      independent <- collapsed$mainPathways
      message(sprintf(
        "collapsePathways: %d sig -> %d independent",
        nrow(sig_dt), length(independent)
      ))
      # Mark non-independent sig pathways as padj = 1 (effectively removes them)
      drop_pw <- setdiff(sig_dt$pathway, independent)
      if (length(drop_pw)) {
        res_dt$padj[res_dt$pathway %in% drop_pw] <- 1
      }
    }

    res <- as.data.frame(res_dt)
    res$database <- classify_database(res$pathway)
    res$contrast <- ctr

    # Jaccard dedup on remaining sig
    sig_after <- res[!is.na(res$padj) & res$padj < padj_cutoff, ]
    sig_dedup <- deduplicate_enrichment(sig_after, pw_list, jaccard_cutoff)
    n_removed <- nrow(sig_after) - nrow(sig_dedup)
    message(sprintf(
      "Jaccard dedup (%.2f): %d -> %d (removed %d)",
      jaccard_cutoff, nrow(sig_after), nrow(sig_dedup), n_removed
    ))

    # Reset padj for terms that didn't survive dedup
    survived <- sig_dedup$pathway
    dedup_drop <- setdiff(sig_after$pathway, survived)
    if (length(dedup_drop)) {
      res$padj[res$pathway %in% dedup_drop] <- 1
    }

    all_results[[ctr]] <- tibble::as_tibble(res)
  }

  long_df <- dplyr::bind_rows(all_results)

  # Union of surviving sig pathways across all contrasts
  sig_union <- unique(long_df$pathway[!is.na(long_df$padj) & long_df$padj < padj_cutoff])
  message(sprintf("\nUnion of sig pathways: %d", length(sig_union)))

  # Filter to union pathways only
  long_df <- long_df[long_df$pathway %in% sig_union, ]

  # Summary
  for (ctr in names(stats_list)) {
    sub <- long_df[long_df$contrast == ctr, ]
    n_sig <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff)
    n_up <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff & sub$NES > 0)
    n_dn <- sum(!is.na(sub$padj) & sub$padj < padj_cutoff & sub$NES < 0)
    message(sprintf("  %s: %d sig (%d up, %d down)", ctr, n_sig, n_up, n_dn))
  }

  list(long_df = long_df, sig_union = sig_union)
}


run_ora_deduplicated <- function(genes, universe, pathways,
                                 em_cutoff = 0.5,
                                 min_size = 10, max_size = 500,
                                 padj_cutoff = 0.05) {
  requireNamespace("fgsea", quietly = TRUE)

  genes <- unique(genes)
  universe <- unique(universe)
  genes <- intersect(genes, universe)

  # Split pathways by database, run fora per database (per-database BH)
  pw_by_db <- split(names(pathways), classify_database(names(pathways)))
  db_results <- list()
  for (db in names(pw_by_db)) {
    db_pw <- pathways[pw_by_db[[db]]]
    if (length(db_pw) < 2) next
    db_res <- fgsea::fora(
      pathways = db_pw,
      genes    = genes,
      universe = universe,
      minSize  = min_size,
      maxSize  = max_size
    )
    db_res <- as.data.frame(db_res)
    if (nrow(db_res) == 0) next
    db_res$database <- db
    db_results[[db]] <- db_res
  }
  res <- do.call(rbind, db_results)
  if (is.null(res) || nrow(res) == 0) {
    return(tibble::tibble(pathway = character()))
  }

  N <- length(universe)
  K <- length(genes)
  res$odds_ratio <- vapply(seq_len(nrow(res)), function(i) {
    a <- res$overlap[i] # hits in pathway
    b <- K - a # foreground not in pathway
    c <- res$size[i] - a # pathway not in foreground
    d <- N - K - c # neither
    if (b == 0 || c == 0) Inf else (a * d) / (b * c)
  }, numeric(1))

  res <- tibble::as_tibble(res)

  # Redundancy is flagged, not dropped, so counts report every significant set
  # while panels draw one representative per cluster. Same split the ranked-list
  # cache uses; coefficient matches build_fgsea_cache.R.
  #
  # Collapsing runs on the FDR < 0.05 subset only. Callers that pass a looser
  # padj_cutoff want the fuller table for plotting, and the similarity loop is
  # quadratic, so flagging thousands of non-significant sets would be both
  # meaningless and slow.
  out <- as.data.frame(res[!is.na(res$padj) & res$padj < padj_cutoff, ])
  if (nrow(out) == 0) {
    out[c("dedup_status", "merged_into", "overlap_jaccard")] <-
      list(character(0), character(0), numeric(0))
    return(tibble::as_tibble(out))
  }
  sig <- out[out$padj < 0.05, ]
  flagged <- em_dedup_report(sig, pathways, em_cutoff)

  out$dedup_status <- NA_character_
  out$merged_into <- NA_character_
  out$overlap_jaccard <- NA_real_
  hit <- match(flagged$pathway, out$pathway)
  out$dedup_status[hit] <- flagged$dedup_status
  out$merged_into[hit] <- flagged$merged_into
  out$overlap_jaccard[hit] <- flagged$overlap_jaccard

  n_kept <- sum(out$dedup_status == "kept", na.rm = TRUE)
  message(sprintf(
    "ORA: %d sig -> %d representatives (%.1f%% redundant)",
    nrow(sig), n_kept,
    if (nrow(sig) > 0) 100 * (nrow(sig) - n_kept) / nrow(sig) else 0
  ))

  tibble::as_tibble(out)
}


classify_database <- function(pathway_names) {
  dplyr::case_when(
    grepl("^HALLMARK_", pathway_names) ~ "Hallmark",
    grepl("^REACTOME_", pathway_names) ~ "Reactome",
    grepl("^KEGG_MEDICUS_", pathway_names) ~ "KEGG",
    grepl("^KEGG_", pathway_names) ~ "KEGG",
    grepl("^GOSLIM_", pathway_names) ~ "GO Slim",
    grepl("^GOBP_", pathway_names) ~ "GO:BP",
    TRUE ~ "Other"
  )
}


# MSigDB pathway ID -> 15 consolidated categories (keyword rules)
CONSOLIDATED_PATHWAY_ORDER <- c(
  "Muscle & Contractile", "Cytoskeleton & Motility", "ECM & Adhesion",
  "Lipid Metabolism", "Carbohydrate & Energy Metabolism",
  "Amino Acid & Cofactor Metabolism",
  "Mitochondria & Energy", "Protein Homeostasis",
  "Transport", "Translation & Ribosome", "Transcription & Chromatin",
  "Immune & Inflammation", "DNA & Cell Cycle", "Circulatory System",
  "Development", "Other"
)

CONSOLIDATED_COLORS <- c(
  "Muscle & Contractile"              = "#E57373",
  "Cytoskeleton & Motility"           = "#FFB74D",
  "ECM & Adhesion"                    = "#FFF176",
  "Lipid Metabolism"                  = "#AED581",
  "Carbohydrate & Energy Metabolism"  = "#81C784",
  "Amino Acid & Cofactor Metabolism"  = "#66BB6A",
  "Mitochondria & Energy"             = "#4DB6AC",
  "Protein Homeostasis"               = "#4FC3F7",
  "Transport"                         = "#7986CB",
  "Translation & Ribosome"            = "#BA68C8",
  "Transcription & Chromatin"         = "#AB47BC",
  "Immune & Inflammation"             = "#A1887F",
  "DNA & Cell Cycle"                  = "#90A4AE",
  "Circulatory System"                = "#CE93D8",
  "Development"                       = "#B0BEC5",
  "Other"                             = "#D0D0D0"
)

classify_pathway_func <- function(ids) {
  rules <- list(
    "Muscle & Contractile"              = "MYOGEN|MYOFIBRIL|SARCOMERE|MUSCLE_|CONTRACTILE|ACTOMYOSIN|MYOSIN|I_BAND",
    "Cytoskeleton & Motility"           = "CYTOSKELET|ACTIN_BIND|STRUCTURAL_MOLECULE|MOTIL|SUPRAMOLECUL",
    "ECM & Adhesion"                    = "EXTRACELLULAR_MATRIX|COLLAGEN|BASEMENT_MEMBRANE|ADHESION|APICAL_JUNCTION|EMT|ENCAPSULATING",
    "Lipid Metabolism"                  = "FATTY_ACID|LIPID|ADIPOGEN|STEROID|SPHINGOLIPID|PHOSPHOLIPID|KETONE",
    "Carbohydrate & Energy Metabolism"  = "GLYCOLY|GLUCONEO|CARBOHYDRATE|PENTOSE|PRECURSOR_METABOL",
    "Amino Acid & Cofactor Metabolism"  = "AMINO_ACID|VITAMIN|COFACTOR|NITROGEN|DETOXIF|DIGEST|XENOBIOT",
    "Mitochondria & Energy"             = "MITOCHOND|OXIDATIVE_PHOSPH|ELECTRON_TRANSFER|RESPIRATORY|OXIDOREDUCT",
    "Protein Homeostasis"               = "PROTEASOM|UBIQUITIN|AUTOPHAGY|MTORC1|PROTEIN_FOLD",
    "Transport"                         = "TRANSPORT(?!.*ELECTRON)|VESICLE|ENDOCYT|SECRETI",
    "Translation & Ribosome"            = "TRANSLAT|RIBOSOM|TRNA|MYC_TARGET",
    "Transcription & Chromatin"         = "TRANSCRIPT|SPLICEOSOM|E2F_TARGET|CHROMATIN|MRNA_PROC",
    "Immune & Inflammation"             = "IMMUN|INFLAMMA|INTERFERON|IL2|IL6|TNFA|NF.KB|COMPLEMENT",
    "DNA & Cell Cycle"                  = "DNA_REPAIR|CELL_CYCLE|MITOTIC|P53_PATHWAY",
    "Circulatory System"                = "ANGIOGEN|BLOOD_VESSEL|HYPOXIA",
    "Development"                       = "UV_RESPONSE|GROWTH_FACTOR|WNT|HEDGEHOG|NOTCH|TGF_BETA|KRAS"
  )
  vapply(toupper(ids), function(id) {
    matches <- character(0)
    for (cat in names(rules)) {
      if (grepl(rules[[cat]], id, perl = TRUE)) matches <- c(matches, cat)
    }
    if (length(matches) > 1) {
      warning(
        "classify_pathway_func: '", id, "' matches multiple categories [",
        paste(matches, collapse = ", "), "]; using first match: ", matches[1]
      )
    }
    if (length(matches) >= 1) {
      return(matches[1])
    }
    if (grepl("METABOL", id, perl = TRUE)) {
      return("Amino Acid & Cofactor Metabolism")
    }
    "Other"
  }, character(1), USE.NAMES = FALSE)
}
