#!/usr/bin/env Rscript
# Build the STRING MCL cluster cache for the nine non-grey WGCNA modules.
# The module labels in Figure 6 have always been described as the module's
# densest STRING cluster, but nothing in the repository ever computed one --
# the labels were read off the ORA output by hand. This script makes the claim
# checkable: it submits each module's gene list to STRING and caches what comes
# back, so the labels can be compared against a real clustering.
#
# Same contract as build_fgsea_cache.R: fetch once, cache to CSV, and let every
# later run read the cache. Once the cache is complete the pipeline needs no
# network at all.

setwd(here::here())

CACHE_PATH <- "04_Figures/F05/c_data/wgcna_string_clusters.csv"
ASSIGN_CSV <- "04_Figures/F05/c_data/wgcna/wgcna_module_assignments.csv"
LABELS_CSV <- "04_Figures/F05/c_data/mod_bio_labels.csv"

# Parameters as reported for the cluster column of S6 Table. One definition,
# used for the fetch and for the provenance stamped into the cache, so the two
# cannot drift apart.
STRING_API <- "https://string-db.org/api/json/network_image_url"
STRING_SPECIES <- 9606L
# MCL inflation controls granularity: higher values cut more edges and return
# more, smaller clusters. The first cache was built at 3.0, which shredded the
# larger modules past the point of meaning -- black's 107 proteins came back as
# 16 clusters whose densest held 16% of the module, and turquoise's 340 as 82.
# That reads as "no coherent biology" when the real answer is that the cut was
# too fine. 1.8 is the coarsest setting that still separates black's two halves
# (a 23-protein translation/UPS cluster and a 17-protein glycolysis cluster)
# instead of merging the module into one blob, so it is the one used here.
STRING_INFLATION <- 1.8
STRING_REQUIRED_SCORE <- 400L
STRING_CALLER <- "YvO_review_wgcna_module_clusters"

# STRING asks callers to stay under one request per second.
STRING_PAUSE_S <- 1

pacman::p_load(dplyr, readr, httr2, digest)

stopifnot(
  "WGCNA module assignments missing - run YvO_WGCNA_run.R" =
    file.exists(ASSIGN_CSV),
  "mod_bio_labels.csv missing - run YvO_WGCNA_run.R" =
    file.exists(LABELS_CSV)
)

labels_df <- read_csv(LABELS_CSV, show_col_types = FALSE)
assign_df <- read_csv(ASSIGN_CSV, show_col_types = FALSE) |>
  filter(module_color != "grey")

module_genes <- assign_df |>
  filter(module_color %in% labels_df$module_color) |>
  distinct(module_color, gene) |>
  arrange(module_color, gene) |>
  (\(d) split(d$gene, d$module_color))()

# Colours are not stable identities across refits, so a cached row is only
# trusted when the module's membership still hashes to what was submitted.
module_hash <- vapply(module_genes, \(g) digest::digest(sort(g)), character(1))

fetch_string_clusters <- function(genes) {
  body <- list(
    identifiers = paste(genes, collapse = "\r"),
    species = STRING_SPECIES,
    required_score = STRING_REQUIRED_SCORE,
    network_clustering_algorithm = "MCL",
    network_clustering_parameter_mcl = STRING_INFLATION,
    caller_identity = STRING_CALLER
  )
  res <- request(STRING_API) |>
    req_body_form(!!!body) |>
    req_retry(max_tries = 3) |>
    req_timeout(300) |>
    req_perform() |>
    # STRING labels its JSON "text/json", which httr2 refuses by default.
    resp_body_json(simplifyVector = TRUE, check_type = FALSE)

  stopifnot("STRING returned no clusters for this module" = NROW(res) > 0)

  tibble(
    cluster_number = as.integer(res$clusterNumber),
    n_proteins = as.integer(res$proteinCount),
    cluster_name = res$clusterNames,
    member_genes = vapply(
      strsplit(res$preferredNames, ",\\s*"),
      \(g) paste(sort(g), collapse = ";"),
      character(1)
    )
  ) |>
    arrange(cluster_number)
}

cached <- if (file.exists(CACHE_PATH)) {
  read_csv(CACHE_PATH, show_col_types = FALSE)
} else {
  NULL
}

fresh_colors <- character(0)
if (!is.null(cached) && nrow(cached)) {
  # Membership is not the only thing a row can go stale against: the clustering
  # parameters are stamped into every row precisely so a change to them can be
  # detected here. Hash alone would have kept the 3.0 clusters forever.
  fresh_colors <- cached |>
    filter(
      module_color %in% names(.env$module_hash),
      gene_hash == unname(.env$module_hash[module_color]),
      inflation == STRING_INFLATION,
      required_score == STRING_REQUIRED_SCORE,
      species == STRING_SPECIES
    ) |>
    pull(module_color) |>
    unique()
}
stale_colors <- setdiff(names(module_genes), fresh_colors)

if (!length(stale_colors)) {
  message(sprintf(
    "STRING cluster cache present (%s) - all %d modules match assignments",
    basename(CACHE_PATH), length(module_genes)
  ))
} else {
  message(sprintf(
    "Fetching STRING MCL clusters (inflation %.1f, score >= %d): %d module(s)",
    STRING_INFLATION, STRING_REQUIRED_SCORE, length(stale_colors)
  ))

  fetched <- lapply(stale_colors, function(mc) {
    genes <- module_genes[[mc]]
    out <- fetch_string_clusters(genes)
    message(sprintf(
      "  %-10s %3d proteins -> %3d clusters (densest %d)",
      mc, length(genes), nrow(out), max(out$n_proteins)
    ))
    Sys.sleep(STRING_PAUSE_S)
    out |>
      mutate(
        module_color = mc,
        module_id = labels_df$module_id[match(mc, labels_df$module_color)],
        n_module_proteins = length(genes),
        fraction_of_module = n_proteins / length(genes),
        gene_hash = unname(module_hash[[mc]]),
        .before = 1
      )
  })

  keep <- if (is.null(cached)) {
    NULL
  } else {
    filter(cached, module_color %in% fresh_colors)
  }

  cache_df <- bind_rows(keep, bind_rows(fetched)) |>
    mutate(
      species = STRING_SPECIES,
      clustering_algorithm = "MCL",
      inflation = STRING_INFLATION,
      required_score = STRING_REQUIRED_SCORE
    ) |>
    arrange(match(module_color, labels_df$module_color), cluster_number)

  cache_df <- cache_df[, c(
    "module_color", "module_id", "n_module_proteins", "cluster_number",
    "n_proteins", "fraction_of_module", "cluster_name", "member_genes",
    "species", "clustering_algorithm", "inflation", "required_score",
    "gene_hash"
  )]

  write_csv(cache_df, CACHE_PATH)
  message(sprintf(
    "Wrote %d cluster rows for %d modules to %s",
    nrow(cache_df), n_distinct(cache_df$module_color), CACHE_PATH
  ))
  cached <- cache_df
}

string_clusters <- read_csv(CACHE_PATH, show_col_types = FALSE)

# Black was clustered by hand on the STRING web interface at this inflation,
# and it is the module the setting was chosen on, so it is the independent
# check that the REST parameters above reproduce what a reader gets from the
# website. Its top two clusters are the whole reason for 1.8: the cache is
# refused if either moves.
black_cluster <- function(clusters, k) {
  is_k <- clusters$module_color == "black" & clusters$cluster_number == k
  row <- clusters[is_k, ]
  if (nrow(row) != 1L) {
    return(list(n = NA_integer_, name = NA_character_))
  }
  list(n = as.integer(row$n_proteins), name = row$cluster_name)
}

matches_web <- function(chk, n, name) {
  identical(chk$n, n) && identical(chk$name, name)
}

n_black_clusters <- sum(string_clusters$module_color == "black")

stopifnot(
  "black no longer resolves into the 11 clusters verified on the STRING site" =
    identical(n_black_clusters, 11L),
  "black cluster 1 no longer matches the verified web result" =
    matches_web(
      black_cluster(string_clusters, 1L),
      23L, "Eukaryotic Translation Elongation"
    ),
  "black cluster 2 no longer matches the verified web result" =
    matches_web(
      black_cluster(string_clusters, 2L),
      17L, "Glycolysis / Gluconeogenesis; Pyruvate metabolism"
    ),
  "the cache holds clusters from a different inflation" =
    all(string_clusters$inflation == STRING_INFLATION),
  "a module colour is missing from the STRING cluster cache" =
    setequal(string_clusters$module_color, names(module_genes))
)

message("STRING cluster cache validated against the verified black clustering")
