# Sourced by F05_data.R after style.R and pathway_utils.R. Writes the hub
# network node table for S6 Table: the Q90-kME hubs of each module, joined by
# their top-decile TOM edges, with kME, gene significance and ORA group.

pacman::p_load(tidyverse, WGCNA, igraph)

allowWGCNAThreads()
set.seed(42)

DAT <- "04_Figures/F05/c_data"

stopifnot(
  "WGCNA gs_phenotype_choices.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/gs_phenotype_choices.csv")),
  "meta.csv missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "meta.csv")),
  "MEs.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "MEs.rds")),
  "kME_all.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "kME_all.rds")),
  "datExpr.rds missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "datExpr.rds")),
  "WGCNA module assignments missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_module_assignments.csv")),
  "WGCNA sft_summary missing — run YvO_WGCNA_run.R first" =
    file.exists(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
)

meta <- read_csv(file.path(DAT, "meta.csv"))
meta$group <- factor(meta$group,
  levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
)
MEs <- readRDS(file.path(DAT, "MEs.rds"))
kME_all <- readRDS(file.path(DAT, "kME_all.rds"))
datExpr <- readRDS(file.path(DAT, "datExpr.rds"))
module_df <- read_csv(file.path(DAT, "wgcna/wgcna_module_assignments.csv"))
sft_csv <- read.csv(file.path(DAT, "wgcna/wgcna_sft_summary.csv"))
NET_POWER <- sft_csv$selected_power[1]

KEY_MODULES <- module_df |>
  filter(module_color != "grey") |>
  count(module_color, sort = TRUE) |>
  pull(module_color)
bg_genes <- unique(module_df$gene)

meta$age_num <- ifelse(meta$age == "Old", 1, 0)

gs_choices <- read_csv(file.path(DAT, "wgcna/gs_phenotype_choices.csv"))
MODULE_GS_PHENO <- setNames(gs_choices$gs_phenotype, gs_choices$module)
gs_label_map <- c(
  age_num = "Age (Young=0, Old=1)",
  delta_VL = "\u0394VL Thickness (cm)",
  delta_LBM = "\u0394Lean Body Mass (kg)",
  VL_thick_cm = "VL Thickness (cm)",
  DXA_LBM_kg = "Lean Body Mass (kg)",
  BMI = "BMI (kg/m\u00b2)"
)
MODULE_GS_LABEL <- setNames(
  ifelse(gs_choices$gs_phenotype %in% names(gs_label_map),
    gs_label_map[gs_choices$gs_phenotype],
    tools::toTitleCase(gsub("_", " ", gs_choices$gs_phenotype))
  ),
  gs_choices$module
)
for (m in KEY_MODULES) {
  if (!(m %in% names(MODULE_GS_PHENO))) {
    MODULE_GS_PHENO[[m]] <- "age_num"
    MODULE_GS_LABEL[[m]] <- "Age (Young=0, Old=1)"
  }
}

if ("delta_VL" %in% MODULE_GS_PHENO && !("delta_VL" %in% colnames(meta))) {
  vl_wide <- meta |>
    filter(!is.na(VL_thick_cm)) |>
    select(subject, time, VL_thick_cm) |>
    pivot_wider(names_from = time, values_from = VL_thick_cm, names_prefix = "VL_") |>
    mutate(delta_VL = VL_Post - VL_Pre)
  meta <- meta |> left_join(vl_wide |> select(subject, delta_VL), by = "subject")
  message(sprintf(
    "  Computed delta_VL for %d subjects (non-NA: %d)",
    nrow(vl_wide), sum(!is.na(vl_wide$delta_VL))
  ))
}

uid2gene <- setNames(module_df$gene, module_df$uniprot_id)

message("Hub protein networks: building the module edge table...")

pw_full <- build_pathway_collection(min_size = 15, include_goslim = FALSE)

select_hubs_q90 <- function(mod) {
  mod_prots <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col <- paste0("kME", mod)
  matched <- intersect(mod_prots, rownames(kME_all))
  mod_kme <- setNames(kME_all[matched, kme_col], matched)
  mod_kme <- mod_kme[!is.na(mod_kme)]
  q90 <- quantile(mod_kme, 0.90)
  names(mod_kme[mod_kme >= q90])
}

assign_groups_ora <- function(gene_names, max_groups = 4, min_group_n = 3) {
  clean_pw_name <- function(name) {
    name |>
      (\(x) gsub("^HALLMARK_|^GOSLIM_|^GOBP_|^REACTOME_|^KEGG_MEDICUS_", "", x))() |>
      (\(x) gsub("_", " ", x))() |>
      str_to_title() |>
      str_trunc(35)
  }

  ora_res <- tryCatch(
    run_ora_deduplicated(genes = gene_names, universe = bg_genes, pathways = pw_full),
    error = function(e) {
      message("  ORA error: ", e$message)
      NULL
    }
  )

  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message("  No significant ORA results")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  ora_res <- ora_res[order(ora_res$padj), ]
  gene_map <- data.frame(gene = character(), pathway = character(), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(ora_res))) {
    hits <- intersect(ora_res$overlapGenes[[i]], gene_names)
    if (length(hits)) {
      gene_map <- rbind(gene_map, data.frame(
        gene = hits, pathway = ora_res$pathway[i],
        stringsAsFactors = FALSE
      ))
    }
  }
  gene_map <- gene_map[!duplicated(gene_map$gene), ]

  term_counts <- table(gene_map$pathway)
  keep <- names(term_counts[term_counts >= min_group_n])
  keep <- head(keep[order(term_counts[keep], decreasing = TRUE)], max_groups)

  if (length(keep) < 2) {
    message("  Fewer than 2 qualifying groups")
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  assignments <- setNames(rep("Other", length(gene_names)), gene_names)
  for (g in gene_names) {
    row <- gene_map[gene_map$gene == g, ]
    if (nrow(row) > 0 && row$pathway[1] %in% keep) {
      assignments[g] <- clean_pw_name(row$pathway[1])
    }
  }
  assignments
}

compute_gs_signed <- function(mod) {
  pheno_col <- MODULE_GS_PHENO[[mod]]
  pheno_vec <- meta[[pheno_col]]
  if (is.null(pheno_vec)) {
    warning(sprintf("  GS phenotype '%s' not found in meta for %s; using age_num", pheno_col, mod))
    pheno_vec <- meta[["age_num"]]
  }
  names(pheno_vec) <- meta$sample_id
  valid_samps <- intersect(meta$sample_id[!is.na(pheno_vec)], rownames(datExpr))
  gs <- cor(datExpr[valid_samps, , drop = FALSE], pheno_vec[valid_samps],
    use = "pairwise.complete.obs"
  )
  setNames(gs[, 1], rownames(gs))
}

hub_node_table <- function(mod) {
  message(toupper(mod))

  hub_ids <- select_hubs_q90(mod)
  hub_genes <- uid2gene[hub_ids]
  hub_genes <- hub_genes[!is.na(hub_genes)]
  hub_ids <- hub_ids[hub_ids %in% names(hub_genes)]
  n_mod <- sum(module_df$module_color == mod)
  message(sprintf("  %d hubs (Q90 of %d)", length(hub_ids), n_mod))

  mod_prots <- intersect(
    module_df$uniprot_id[module_df$module_color == mod],
    colnames(datExpr)
  )
  adj_mod <- adjacency(datExpr[, mod_prots], power = NET_POWER, type = "signed hybrid")
  tom_mod <- TOMsimilarity(adj_mod, TOMType = "signed")
  colnames(tom_mod) <- rownames(tom_mod) <- mod_prots

  tom_sub <- tom_mod[hub_ids, hub_ids]
  tom_q90 <- quantile(tom_sub[upper.tri(tom_sub)], 0.90)

  g <- graph_from_adjacency_matrix(tom_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- delete_edges(g, which(E(g)$weight < tom_q90))
  iso <- which(degree(g) == 0)
  if (length(iso)) g <- delete_vertices(g, iso)

  node_uids <- V(g)$name
  node_genes <- uid2gene[node_uids]
  node_gs <- compute_gs_signed(mod)[node_uids]
  node_gs[is.na(node_gs)] <- 0
  groups <- assign_groups_ora(node_genes)

  data.frame(
    name = node_uids, gene = node_genes,
    kME = setNames(kME_all[node_uids, paste0("kME", mod)], node_uids),
    GS = node_gs, func_grp = groups[node_genes]
  )
}

node_data_list <- setNames(lapply(KEY_MODULES, hub_node_table), KEY_MODULES)
all_node_df <- bind_rows(lapply(KEY_MODULES, function(mod) {
  nd <- node_data_list[[mod]]
  if (is.null(nd) || nrow(nd) == 0) {
    return(NULL)
  }
  tibble(
    uniprot_id = nd$name, module = mod, gene = nd$gene,
    kME = nd$kME, GS = nd$GS, functional_group = nd$func_grp,
    GS_pheno = MODULE_GS_LABEL[[mod]]
  )
}))
write_csv(all_node_df, file.path(DAT, "04_panel_D_hub_network.csv"))

message(sprintf("  Hub edge table written for %d modules", length(KEY_MODULES)))
