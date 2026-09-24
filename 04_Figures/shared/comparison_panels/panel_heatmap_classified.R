# Shared classified-heatmap engine: row-z-scored abundance, rows grouped
# directly by GO Slim + Hallmark functional classification (membership-based,
# not enrichment or expression clustering: every protein's primary category
# is its smallest matching set across both databases pooled, which maximizes
# how many proteins get a real label instead of falling to "Unclassified").
# Categories under `min_group_size` fold into "Other classified"; genuinely
# unclassified proteins form their own group, always last. Within a group,
# rows sort by the primary log2FC column (descending) -- there is no
# expression-pattern clustering step. Gene names render bold when classified,
# plain grey when not, so the labeling coverage is visible on the plot
# itself, not just claimed in a caption. Config-driven, following
# panel_D_nes_scatter.R's calling convention: set `cfg`,
# then source this file. Produces `pH_heat` (ggdraw composite), `pH_zmat`,
# `pH_group_id`, `pH_K` (group count), `pH_cluster_membership`.
#
# cfg fields:
#   protein_df       pre-filtered data frame (one row per displayed protein),
#                     must carry `gene`, every column in `sample_cols`, and
#                     every column named in `logfc_cols`/`pi_col`
#   sample_cols       character vector of Pre/Post sample columns
#   logfc_cols        named character vector; names are display labels, first
#                      entry drives within-group row order and legend title
#   pi_col            column name for a Pi < 0.05 marker annotation, or NULL
#   fdr_col           column name for an adj.P.Val < 0.05 marker annotation,
#                      or NULL -- symmetric to pi_col, for panels gated on Pi
#                      that still want to show which rows also clear FDR
#   min_group_size    categories with fewer members fold into "Other
#                      classified" (default 4)
#   cluster_within_group  if TRUE, rows within each functional group are
#                      ordered by real hierarchical clustering (Euclidean,
#                      Ward.D2, on the row-z-scored matrix) instead of by
#                      primary log2FC (default FALSE)
#   universe_genes    character vector, full gene universe (kept for parity
#                      with callers; not used for classification itself)
#   title, subtitle_fmt (sprintf(subtitle_fmt, n_total, K))
#   file_stub         output filename stem, e.g. "aging_heatmap"
#   rpt_png, rpt_pdf, dat_ht   output directories (already created)
#   comp_w, comp_h    composite panel size in mm

pacman::p_load(
  dplyr, tidyr, tibble, stringr, purrr, ggplot2, cowplot,
  ComplexHeatmap, circlize, grid
)

pH_pretty_term <- function(x) {
  x |>
    str_remove("^(GOBP|GOCC|GOMF|GOSLIM|REACTOME|KEGG(_MEDICUS)?(_REFERENCE)?|HALLMARK|WP)_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_replace("complex i$", "complex I") |>
    str_replace("^Trna", "tRNA") |>
    str_replace("^Rrna", "rRNA") |>
    str_replace("^Mrna", "mRNA")
}

pH_df <- cfg$protein_df
pH_min_group_size <- cfg$min_group_size %||% 4

pH_abund <- as.matrix(pH_df[, cfg$sample_cols])
rownames(pH_abund) <- pH_df$gene
pH_zmat <- t(scale(t(pH_abund)))

pH_sample_block <- factor(
  paste(
    if_else(grepl("^Y", cfg$sample_cols), "Young", "Old"),
    sub(".*_", "", cfg$sample_cols)
  ),
  levels = c("Young Pre", "Young Post", "Old Pre", "Old Post")
)

# Functional classification, not enrichment: GO Slim + Hallmark are curated,
# non-redundant vocabularies built for exactly this kind of readable
# categorization (that's what "Slim" means), pooled across both databases so
# a protein missing from one still gets classified via the other -- no
# p-value or background comparison involved.
pH_class_pool <- build_pathway_collection(
  min_size = 5, max_size = 800,
  include_goslim = TRUE, exclude_variants = TRUE
)
pH_class_pool <- pH_class_pool[classify_database(names(pH_class_pool)) %in% c("Hallmark", "GO Slim")]

pH_membership <- tibble(
  pathway = rep(names(pH_class_pool), lengths(pH_class_pool)),
  gene = unlist(pH_class_pool, use.names = FALSE)
) |>
  filter(gene %in% pH_df$gene)

# Primary category per protein = whichever matching term captures the most
# *other* proteins in this specific set, not globally-smallest-set-size.
# Picking by global specificity scatters proteins across the long tail of
# narrow terms, half of them singletons; picking by local size concentrates
# them into fewer, more populated groups covering roughly half again as many
# proteins. The point of grouping here is coverage for display, not maximal
# per-protein specificity. Proteins matching nothing in the pool
# stay "Unclassified" rather than being silently dropped from the count.
pH_local_n <- pH_membership |> count(pathway, name = "local_n")
pH_primary <- pH_membership |>
  left_join(pH_local_n, by = "pathway") |>
  slice_max(local_n, by = gene, n = 1, with_ties = FALSE) |>
  select(gene, pathway)

pH_lfc_primary <- pH_df[[cfg$logfc_cols[1]]]
names(pH_lfc_primary) <- pH_df$gene

pH_cluster_membership <- tibble(gene = pH_df$gene) |>
  left_join(pH_primary, by = "gene") |>
  mutate(pathway = coalesce(pathway, "Unclassified"))

# Categories below min_group_size fold into "Other classified" -- keeps the
# heatmap from fragmenting into dozens of one- or two-protein blocks while
# every classified protein still renders bold and is still named in the
# exported CSV.
pH_cat_counts <- pH_cluster_membership |>
  filter(pathway != "Unclassified") |>
  count(pathway, sort = TRUE)
pH_big_cats <- pH_cat_counts$pathway[pH_cat_counts$n >= pH_min_group_size]

pH_cluster_membership <- pH_cluster_membership |>
  mutate(
    group_key = case_when(
      pathway == "Unclassified" ~ "Unclassified",
      pathway %in% pH_big_cats ~ pathway,
      TRUE ~ "Other classified"
    ),
    classified = pathway != "Unclassified"
  )

# Group order: real categories largest-first, "Other classified" next,
# "Unclassified" always last regardless of size -- signal before residuals.
pH_group_order <- c(
  pH_cat_counts$pathway[pH_cat_counts$pathway %in% pH_big_cats],
  "Other classified", "Unclassified"
)
pH_group_order <- pH_group_order[
  pH_group_order %in% pH_cluster_membership$group_key
]
pH_K <- length(pH_group_order)

pH_group_n <- table(pH_cluster_membership$group_key)
# Row-title height is set by the slice's row count, not by label length --
# a 4-line wrapped name on a 3-row slice overflows into its neighbors. Cap
# every label at 2 lines (name, then "(n = X)") and truncate the name if it
# still doesn't fit, so small groups with long term names stay legible.
pH_group_lab <- setNames(
  vapply(pH_group_order, function(g) {
    name <- if (g %in% c("Other classified", "Unclassified")) g else pH_pretty_term(g)
    if (nchar(name) > 24) name <- paste0(substr(name, 1, 21), "...")
    sprintf("%s\n(n = %d)", name, pH_group_n[[g]])
  }, character(1)),
  pH_group_order
)
pH_group_confident <- setNames(
  !pH_group_order %in% c("Other classified", "Unclassified"),
  pH_group_order
)

pH_cluster_membership <- pH_cluster_membership |>
  mutate(
    group_label = pH_group_lab[group_key],
    lfc = pH_lfc_primary[gene]
  ) |>
  arrange(factor(group_key, levels = pH_group_order), desc(lfc))

# Reorder the display matrix to match: grouped by category (largest first),
# sorted by primary log2FC within each group -- no expression clustering.
pH_zmat <- pH_zmat[pH_cluster_membership$gene, , drop = FALSE]
pH_df <- pH_df[match(pH_cluster_membership$gene, pH_df$gene), ]

pH_row_split <- factor(pH_cluster_membership$group_label, levels = unname(pH_group_lab[pH_group_order]))
pH_title_gp <- gpar(
  fontsize = 6.5,
  fontface = ifelse(pH_group_confident[pH_group_order], "bold", "italic"),
  col = ifelse(pH_group_confident[pH_group_order], "black", "grey45")
)

pH_lfc_cap <- max(abs(unlist(pH_df[, unname(cfg$logfc_cols)])), na.rm = TRUE)
pH_lfc_fun <- colorRamp2(c(-pH_lfc_cap, 0, pH_lfc_cap), c("#4393C3", "white", "#D6604D"))

pH_row_anno_args <- list()
for (nm in names(cfg$logfc_cols)) {
  pH_row_anno_args[[nm]] <- pH_df[[cfg$logfc_cols[[nm]]]]
}
pH_col_args <- setNames(rep(list(pH_lfc_fun), length(cfg$logfc_cols)), names(cfg$logfc_cols))
pH_show_legend <- c(TRUE, rep(FALSE, length(cfg$logfc_cols) - 1))
pH_legend_param <- setNames(
  list(list(
    title = "log2FC", direction = "horizontal",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5)
  )),
  names(cfg$logfc_cols)[1]
)

pH_pi_col <- c(`Π < 0.05` = "grey20", `n.s.` = "grey90")
if (!is.null(cfg$pi_col)) {
  pH_row_anno_args[["Π sig."]] <- if_else(pH_df[[cfg$pi_col]] < 0.05, "Π < 0.05", "n.s.")
  pH_col_args[["Π sig."]] <- pH_pi_col
  pH_show_legend <- c(pH_show_legend, TRUE)
  pH_legend_param[["Π sig."]] <- list(
    title = "Π < 0.05", direction = "horizontal",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5)
  )
}

pH_fdr_col <- c(`FDR < 0.05` = "grey20", `n.s.` = "grey90")
if (!is.null(cfg$fdr_col)) {
  pH_row_anno_args[["FDR sig."]] <- if_else(pH_df[[cfg$fdr_col]] < 0.05, "FDR < 0.05", "n.s.")
  pH_col_args[["FDR sig."]] <- pH_fdr_col
  pH_show_legend <- c(pH_show_legend, TRUE)
  pH_legend_param[["FDR sig."]] <- list(
    title = "FDR < 0.05", direction = "horizontal",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5)
  )
}

ha_row <- do.call(rowAnnotation, c(
  pH_row_anno_args,
  list(
    col = pH_col_args,
    annotation_name_gp = gpar(fontsize = 5),
    annotation_name_side = "top",
    annotation_name_rot = 45,
    width = unit(9, "mm"), gap = unit(0.6, "mm"),
    show_legend = pH_show_legend,
    annotation_legend_param = pH_legend_param
  )
))

pH_z_cap <- 2.5
pH_z_fun <- colorRamp2(c(-pH_z_cap, 0, pH_z_cap), c("#4393C3", "white", "#D6604D"))
# The cap keeps a dense composite panel from turning into a wall of gene names.
# A standalone supplementary list is the opposite case -- the names are the
# point -- so both the cap and the label size are callable.
pH_name_cap <- if (is.null(cfg$row_name_cap)) 200 else cfg$row_name_cap
pH_name_pt <- if (is.null(cfg$row_name_pt)) 3.8 else cfg$row_name_pt
pH_show_names <- nrow(pH_zmat) <= pH_name_cap
pH_row_names_gp <- gpar(
  fontsize = pH_name_pt,
  fontface = ifelse(pH_cluster_membership$classified, "bold.italic", "italic"),
  col = ifelse(pH_cluster_membership$classified, "black", "grey50")
)
pH_cluster_within <- isTRUE(cfg$cluster_within_group)

# NA-safe distance: some proteins have no observed value for a given sample
# (dropout survives imputation upstream in rare cases), which would otherwise
# make dist() return NA for that pair. Clustering fills with the protein's
# row mean; the displayed cells keep their true NA (grey92) regardless.
pH_fill_dist <- function(m) {
  i <- which(is.na(m))
  if (length(i)) m[i] <- rowMeans(m, na.rm = TRUE)[row(m)[i]]
  dist(m)
}

ht <- Heatmap(
  pH_zmat,
  name = "z", col = pH_z_fun, na_col = "grey92",
  cluster_rows = pH_cluster_within, cluster_columns = FALSE,
  clustering_distance_rows = pH_fill_dist, clustering_method_rows = "ward.D2",
  show_row_dend = pH_cluster_within, row_dend_width = unit(6, "mm"),
  show_row_names = pH_show_names,
  row_names_side = "left",
  row_names_gp = pH_row_names_gp,
  row_split = pH_row_split, column_split = pH_sample_block,
  cluster_row_slices = FALSE,
  row_gap = unit(1.5, "mm"), column_gap = unit(3.2, "mm"),
  show_column_names = FALSE,
  column_title_gp = gpar(fontsize = 7, fontface = "bold"),
  row_title_gp = pH_title_gp,
  row_title_rot = 0,
  right_annotation = ha_row,
  heatmap_legend_param = list(
    title = "abundance (row z-score)", direction = "horizontal",
    title_gp = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5)
  )
)

pH_draw_ht <- function() {
  draw(ht,
    heatmap_legend_side = "bottom", merge_legend = TRUE,
    padding = unit(c(1, 2, 8, 12), "mm")
  )
}
pH_g_ht <- grid.grabExpr(pH_draw_ht(), wrap = TRUE, wrap.grobs = TRUE)

# The base png and pdf devices route text through mbcsToSbcs, which drops
# every Pi in these labels: the panel PDFs shipped reading " sig." rather than
# "Pi sig.". agg_png and open_pdf keep the glyph.
ragg::agg_png(
  file.path(cfg$rpt_png, sprintf("MAIN_panel_%s.png", cfg$file_stub)),
  width = cfg$comp_w, height = cfg$comp_h - 12, units = "mm", res = 300,
  background = "white"
)
pH_draw_ht()
dev.off()
open_pdf(
  file.path(cfg$rpt_pdf, sprintf("MAIN_panel_%s.pdf", cfg$file_stub)),
  cfg$comp_w / 25.4, (cfg$comp_h - 12) / 25.4
)
pH_draw_ht()
dev.off()

pH_txt <- composite_text_sizes(cfg$comp_w)
pH_subtitle_str <- str_wrap(sprintf(cfg$subtitle_fmt, nrow(pH_zmat), pH_K), width = 110)
pH_sub_lines <- lengths(gregexpr("\n", pH_subtitle_str)) + 1
pH_header_h <- 0.958 - 0.014 * max(0, pH_sub_lines - 1)
pH_heat <- ggdraw() +
  draw_grob(pH_g_ht, x = 0, y = 0, width = 1, height = pH_header_h) +
  draw_label(cfg$title,
    x = 0.03, y = 0.99, size = pH_txt$title, fontface = "bold",
    hjust = 0, vjust = 1
  ) +
  draw_label(
    pH_subtitle_str,
    x = 0.03, y = 0.972, size = pH_txt$subtitle, fontface = "bold.italic",
    hjust = 0, vjust = 1, colour = "grey40", lineheight = 0.9
  )

heat_df <- as.data.frame(pH_zmat) |>
  mutate(gene = rownames(pH_zmat), group = pH_cluster_membership$group_key, .before = 1)
write_csv(heat_df, file.path(cfg$dat_ht, sprintf("%s_groups.csv", cfg$file_stub)))
write_csv(
  pH_cluster_membership,
  file.path(cfg$dat_ht, sprintf("%s_cluster_classification.csv", cfg$file_stub))
)

pH_n_classified <- sum(pH_cluster_membership$classified)
message(sprintf(
  "%s heatmap: %d proteins x %d samples, %d functional groups, %d/%d classified (%s)",
  cfg$fig_id, nrow(pH_zmat), ncol(pH_zmat), pH_K,
  pH_n_classified, nrow(pH_zmat), paste(pH_group_n[pH_group_order], collapse = "/")
))
