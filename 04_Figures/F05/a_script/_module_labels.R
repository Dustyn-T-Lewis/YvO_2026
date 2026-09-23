# Two sheets behind the module names on Figure 6 panel A: the protein family
# that dominates each module's kME core, and the GO:BP ancestor its enriched
# terms roll up to.

pacman::p_load(dplyr, readr, stringr)
# Loaded but not attached: AnnotationDbi exports select(), and masking dplyr's
# breaks every downstream sheet that builds a data frame.
stopifnot(requireNamespace("GO.db", quietly = TRUE))

# kME >= 0.6 is the module's own definition of a core: a protein correlating
# that strongly with its eigengene is driving the module rather than riding it.
#
# 0.7 was too tight to label the small modules. It left black with 11 core
# proteins against green's 58, so the same rule gave one module five times the
# evidence of another and the thin ones produced the shakiest names. Two labels
# looked unevidenced purely because of it: black's glycolysis half has no
# glycolytic enzyme above 0.7 and gains PFKM, PFKL and STBD1 below it, and
# red's 60S call rests on 28 RPL proteins at 0.7 and 36 at 0.6. No module's
# dominant family changes between the two cuts; every one is reinforced.
KME_CORE <- 0.6

# The ancestor must cover a fifth of the module's significant terms and at least
# three of them. Both floors were read off the descent, not chosen first:
#   - at 0.25 brown stops at "metabolic process" (depth 2), the useless band;
#   - at 0.15 blue follows its artefact branch down to "DNA metabolic process",
#     which is telomerase and reverse-transcriptase terms in a chaperone module.
# 0.20 is the only value between the two where no module lands at depth <= 2 and
# none lands on a branch its core family contradicts. The count floor matters
# because black and magenta have ten significant terms each, so one term is 10%
# of the module: without it black descends to "glycolytic process through
# fructose-6-phosphate" at depth 13 on the strength of two terms.
ORA_COVER_FRAC <- 0.20
ORA_COVER_MIN <- 3L

BP_ROOT <- "GO:0008150"
# A "regulation of X" head names a relationship to a process, not a process, so
# it cannot fill the label's second slot. Barred as a head only -- regulation
# terms stay in the denominator, where they count against any head that claims
# the module.
BP_REGULATION <- "GO:0065007"

.go_bp_term <- local({
  terms <- AnnotationDbi::Term(GO.db::GOTERM)
  terms[AnnotationDbi::Ontology(GO.db::GOTERM) == "BP"]
})

.go_ancestor <- as.list(GO.db::GOBPANCESTOR)
.go_child <- as.list(GO.db::GOBPCHILDREN)
.go_parent <- as.list(GO.db::GOBPPARENTS)

# MSigDB derives its GOBP_ names from the GO term name, so normalising both
# sides the same way recovers the GO id without an extra mapping file.
# str_to_title() lowercases the acronyms GO spells out in full caps, and "Atp
# Synthesis" on a figure reads as a typo.
title_case_term <- function(x) {
  str_replace_all(
    str_to_title(x),
    c(
      "\\bAtp\\b" = "ATP", "\\bAdp\\b" = "ADP", "\\bGtp\\b" = "GTP",
      "\\bDna\\b" = "DNA", "\\bRna\\b" = "RNA", "\\bMrna\\b" = "mRNA",
      "\\bRrna\\b" = "rRNA", "\\bNadh\\b" = "NADH", "\\bNad\\b" = "NAD"
    )
  )
}

msigdb_key <- function(x) {
  paste0("GOBP_", str_remove_all(
    str_replace_all(str_to_upper(x), "[^A-Za-z0-9]+", "_"), "^_|_$"
  ))
}

.msigdb_to_go <- local({
  key <- msigdb_key(.go_bp_term)
  out <- setNames(names(.go_bp_term), key)
  out[!duplicated(names(out))]
})

# Longest path to the root, so a term that is reachable both ways is scored by
# the more specific of its two positions.
.depth_cache <- new.env(parent = emptyenv())
go_depth <- function(id) {
  if (identical(id, BP_ROOT)) {
    return(0L)
  }
  if (!is.null(.depth_cache[[id]])) {
    return(.depth_cache[[id]])
  }
  parents <- setdiff(.go_parent[[id]], c("all", NA))
  # Seeded before the recursion so a cycle in the graph cannot loop forever.
  .depth_cache[[id]] <- 1L
  d <- if (length(parents) == 0L) {
    1L
  } else {
    1L + max(vapply(parents, go_depth, integer(1)))
  }
  .depth_cache[[id]] <- d
  d
}

has_ancestor <- function(ids, target) {
  vapply(ids, function(x) target %in% .go_ancestor[[x]], logical(1))
}

go_depth_of <- function(ids) vapply(ids, go_depth, integer(1))

is_regulation <- function(ids) {
  ids == BP_REGULATION | has_ancestor(ids, BP_REGULATION)
}

# Every ancestor of every significant term, scored by how much of the module it
# accounts for.
ancestor_coverage <- function(go_ids, padj) {
  cand <- unique(unlist(lapply(
    go_ids, function(g) c(g, setdiff(.go_ancestor[[g]], "all"))
  )))
  cand <- cand[cand != BP_ROOT]
  bind_rows(lapply(cand, function(a) {
    hit <- go_ids == a | has_ancestor(go_ids, a)
    tibble(
      go_id = a,
      n_covered = sum(hit),
      best_padj = min(padj[hit])
    )
  })) |>
    mutate(
      depth = go_depth_of(go_id),
      term = as.character(.go_bp_term[go_id])
    )
}

# Walk down from the root along the child that keeps the most of the module,
# and stop at the step where that coverage collapses. Taking the deepest term
# clearing the same floor instead lets a minor branch with many redundant terms
# win the label; descending keeps the head on the module's dominant path.
descend_to_ancestor <- function(cover, n_terms) {
  open <- cover |> filter(!is_regulation(go_id))
  node <- BP_ROOT
  repeat {
    kids <- intersect(.go_child[[node]], open$go_id)
    if (length(kids) == 0L) break
    best <- open |>
      filter(go_id %in% kids) |>
      arrange(desc(n_covered), desc(depth), best_padj, go_id) |>
      slice_head(n = 1)
    collapsed <- best$n_covered < ORA_COVER_MIN ||
      best$n_covered / n_terms < ORA_COVER_FRAC
    if (collapsed) break
    node <- best$go_id
  }
  if (identical(node, BP_ROOT)) NULL else filter(open, go_id == node)
}

.mod_labels <- read_csv(f06("mod_bio_labels.csv"), show_col_types = FALSE)
.kme <- read_csv(f06("wgcna_kme_all.csv"), show_col_types = FALSE)

# The relaxed background, but the choice is free: at padj < 0.05 the strict and
# relaxed GO:BP sets hold the same pathways with the same adjusted p-values in
# all nine modules. Relaxed is the superset of terms *tested*, and is what the
# triptych panels already draw on.
.enrich_relaxed <- read_csv(
  f06("03_panel_B_triptych_enrichment.csv"),
  show_col_types = FALSE
)

# Largest protein family in each module's core, read off the KME_CORE lists on
# 2026-09-16. Curated rather than matched by symbol prefix: the ribosomal,
# proteasomal and respiratory families fall out of a prefix cleanly but the
# sarcomere and beta-oxidation ones do not, and a rule tuned until it agreed
# would be the curation with extra steps.
core_family_lookup <- c(
  turquoise = "Beta-Oxidation", blue = "Chaperones", brown = "20S Proteasome",
  yellow = "Sarcomere", green = "Respiratory Chain",
  black = "Translation Factors",
  red = "60S Ribosome", pink = "ATP Synthase", magenta = "40S Ribosome"
)
# Stem of each core family, for the tie-break below: it is how "20S Proteasome"
# recognises "Proteasomal Ubiquitin Independent Protein Catabolic Process" as
# its own biology. Curated beside the family names rather than clipped off them,
# because "Chaperones" and "Respiratory Chain" have no stem a rule would find.
core_family_stem <- c(
  turquoise = "Oxidation", blue = "Chaperon", brown = "Proteasom",
  yellow = "Sarcomere", green = "Respirat", black = "Translation Factor",
  red = "Ribosom", pink = "ATP Synth", magenta = "Ribosom"
)
core_family_n <- c(
  turquoise = 15L, blue = 10L, brown = 6L, yellow = 27L, green = 48L,
  black = 4L, red = 28L, pink = 9L, magenta = 14L
)

.module_core <- .kme |>
  filter(kME >= KME_CORE) |>
  arrange(module_color, desc(kME)) |>
  summarise(
    core_n = n(),
    top_core_genes = paste(head(gene, 10), collapse = ", "),
    .by = module_color
  ) |>
  right_join(
    .mod_labels[, c("module_color", "module_id", "n_proteins")],
    by = "module_color"
  ) |>
  mutate(
    core_family = unname(core_family_lookup[module_color]),
    family_n = unname(core_family_n[module_color]),
    family_pct = round(100 * family_n / core_n)
  ) |>
  arrange(module_id) |>
  dplyr::select(
    module_color, module_id, n_proteins, core_n, core_family, family_n,
    family_pct, top_core_genes
  ) |>
  as.data.frame()

stopifnot(
  "a module has no curated core family - re-read the kME core lists" =
    !anyNA(.module_core$core_family),
  "a curated core family is larger than the core it came from" =
    all(.module_core$family_n <= .module_core$core_n)
)

# THE RULE: each module's label takes its most significant GO:BP term strictly
# beneath its roll-up ancestor, skipping any term another module has already
# been given, and falling back to the module's best term overall when the branch
# has nothing left to give.
#
# The ancestor is what keeps the label off the artifacts -- it picks the branch
# the module's enrichment actually sits on, so red never reaches "Negative
# Regulation Of Syncytium Formation" (padj 9.6e-41, overlap entirely ribosomal)
# and blue never reaches the GTSE1 term. Taking the ancestor itself as the label
# was tried first and climbed too far: blue landed on "Macromolecule Metabolic
# Process", brown one step from the root. The descendant carries the same
# branch at a depth a reader can use.
#
# Green skips its own branch. That branch is ATP and nucleotide synthesis, and
# the ATP terms are pink's half of the process: green builds the proton gradient
# with 48 of its 58 core proteins from complexes I, III and IV, pink spends it
# with 9 of 22 from the ATP synthase. Green scores better on ATP synthesis than
# pink does, so no p-value ordering splits them -- it is stated here instead.
VETO_OWN_BRANCH <- "green"

# Where the leading candidates are separated by less than this ratio in adjusted
# p-value, the one whose biology matches the module's core family wins. The
# ordering inside such a margin is noise and the tie-break says so out loud
# instead of letting the noise name the module. Brown is the case it exists for:
# its best term is Monosaccharide Biosynthetic Process at padj 0.0097 and
# Proteasomal Ubiquitin Independent Protein Catabolic Process is third at 0.010,
# a gap of 0.0003, and picking the first left a label whose two halves pointed
# opposite ways. Brown's whole profile is weak -- six more terms sit inside
# padj 0.023, among them Cilium Movement and Cilium Or Flagellum Dependent Cell
# Motility, neither of which belongs in skeletal muscle -- which is the reason
# to distrust the ordering rather than the terms.
ORA_TIE_RATIO <- 1.1

# One label is set by hand. Turquoise's rule answer is Alpha Amino Acid
# Metabolic Process at padj 9e-10, because its branch is carboxylic-acid
# catabolism and the fatty-acid half of that is what "Beta-Oxidation" already
# says. Naming the amino-acid half instead gives a label with two fuels and no
# pathway. The TCA cycle is where both arrive, and it is the module's 18th term
# by p-value and outside the branch, so no coverage rule reaches it.
LABEL_OVERRIDE <- c(turquoise = "Tricarboxylic Acid Cycle")

.module_ora <- local({
  sig <- .enrich_relaxed |>
    filter(database == "GO:BP", padj < 0.05) |>
    mutate(go_id = unname(.msigdb_to_go[pathway])) |>
    filter(!is.na(go_id))

  # Black returns nothing at padj < 0.05. That is the module's result and not a
  # gap to paper over -- an 11-protein core that fails to reform in most LOSO
  # refits has little to annotate -- so its ten nearest misses stand in and the
  # significant column says FALSE.
  module_terms <- function(m) {
    s <- arrange(filter(sig, module == m), padj, desc(go_depth_of(go_id)))
    if (nrow(s) > 0L) {
      return(s)
    }
    .enrich_relaxed |>
      filter(module == m, database == "GO:BP") |>
      mutate(go_id = unname(.msigdb_to_go[pathway])) |>
      filter(!is.na(go_id)) |>
      arrange(padj, desc(go_depth_of(go_id))) |>
      slice_head(n = 10)
  }

  prepared <- lapply(.mod_labels$module_color, function(m) {
    s <- module_terms(m)
    anc_head <- descend_to_ancestor(ancestor_coverage(s$go_id, s$padj), nrow(s))
    under <- filter(
      s, go_id != anc_head$go_id, has_ancestor(go_id, anc_head$go_id)
    )
    pool <- if (m %in% VETO_OWN_BRANCH) s[0, ] else under
    list(module = m, terms = s, anc = anc_head, under = under, pool = pool)
  })
  names(prepared) <- .mod_labels$module_color

  # A shared term goes to the module that scores better on it, so the modules
  # are resolved best-p-value first and each claim is final.
  best_p <- vapply(prepared, function(x) {
    min(c(x$pool$padj, x$terms$padj))
  }, numeric(1))

  taken <- character(0)
  chosen <- list()
  for (m in names(sort(best_p))) {
    x <- prepared[[m]]
    forced <- LABEL_OVERRIDE[m]
    pick <- if (!is.na(forced)) {
      filter(x$terms, Description == unname(forced))
    } else {
      free <- filter(x$pool, !Description %in% taken)
      if (nrow(free) == 0L) free <- filter(x$terms, !Description %in% taken)
      tied <- filter(free, padj <= min(padj) * ORA_TIE_RATIO)
      own <- filter(tied, str_detect(Description, fixed(core_family_stem[[m]])))
      slice_head(if (nrow(own) > 0L) own else free, n = 1)
    }
    taken <- c(taken, pick$Description)
    chosen[[m]] <- mutate(
      pick[, c("Description", "go_id", "padj")],
      module_color = m,
      source = if (!is.na(forced)) {
        "set by hand"
      } else if (pick$go_id %in% x$under$go_id) {
        paste0("descendant of ", title_case_term(x$anc$term))
      } else {
        paste0("outside ", title_case_term(x$anc$term))
      },
      n_terms_rolled_up = paste0(x$anc$n_covered, "/", nrow(x$terms)),
      top_descendant_terms = paste(
        head(x$under$Description, 3),
        collapse = "; "
      )
    )
  }

  bind_rows(chosen) |>
    left_join(
      .mod_labels[, c("module_color", "module_id")],
      by = "module_color"
    ) |>
    arrange(module_id) |>
    transmute(
      module_color, module_id,
      label_term = title_case_term(Description),
      label_go_id = go_id,
      label_depth = go_depth_of(go_id),
      ancestor_or_descendant = source,
      n_terms_rolled_up,
      best_padj = padj,
      significant = padj < 0.05,
      top_descendant_terms
    ) |>
    as.data.frame()
})

stopifnot(
  "a module has no label term" = nrow(.module_ora) == nrow(.mod_labels),
  "two modules were given the same GO:BP term" =
    !anyDuplicated(.module_ora$label_term),
  "a label term came back too general to print" =
    all(.module_ora$label_depth >= 5)
)
