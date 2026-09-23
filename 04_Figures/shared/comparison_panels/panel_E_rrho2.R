# Shared RRHO2 Panel E — used by F04 (concordance) and F05 (reversal)

source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_apply.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, ggrepel, patchwork, msigdbr, fgsea, RRHO2)

PE_W <- cfg$panel_w %||% 89
# Descriptive mode draws the overlap map with hotspot sizes only: no corner
# test and no quadrant ORA. It is for contrast pairs that share a group mean
# (Aging and Training(Old) share Old_Pre), where any p-value on the overlap
# would be inflated by the shared noise.
descriptive <- isTRUE(cfg$descriptive)
file_stem <- cfg$file_stem %||% "MAIN_panel_E_rrho2"

RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

# Supplementary output dirs (F05 only)
if (!is.null(cfg$supp)) {
  dir.create(cfg$supp$rpt_png, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfg$supp$rpt_pdf, recursive = TRUE, showWarnings = FALSE)
}

pdf_device <- get_pdf_device()

dep_df <- read_csv(DEP_RESULTS, show_col_types = FALSE)

rr_df <- dep_df |>
  transmute(gene,
    t_1 = .data[[cfg$t_col_1]],
    t_2 = .data[[cfg$t_col_2]]
  ) |>
  filter(!is.na(t_1) & !is.na(t_2)) |>
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# One hypergeometric test per quadrant at a single prespecified depth.
#
# The RRHO heatmap's summit is a maximum over the whole rank grid, so it has no
# hypergeometric null and cannot be quoted as a p-value (Plaisier et al. 2010,
# NAR 38:e169; Cahill et al. 2018, Sci Rep 8:9588). Fixing the depth in advance
# restores an exact test that needs no multiplicity correction. A fixed
# FDR-threshold overlap is not an option here: Training_Old returns zero
# FDR-significant proteins in the limma tree, which collapses the 2x2 table
# (Irizarry et al. 2005, Nat Methods 2:345, for the top-k concordance idiom).
# A fraction rather than a count: the depth is then one stated rule that ports
# to any pair of lists, instead of a number that has to be re-justified when n
# changes. 10% is a round convention, not a prescription from the RRHO papers,
# which report no corner test at all. Shallower fails: at 5% the discordant
# corners hold zero proteins against 5.3 expected and the depletion falls to
# p = 0.03.
TOP_FRAC <- 0.10
TOP_K <- round(TOP_FRAC * n_shared)

quadrant_hyper <- function(dir_1, dir_2, k = TOP_K) {
  pick <- function(score, dir) {
    which(rank(if (dir > 0) -score else score, ties.method = "first") <= k)
  }
  a <- pick(rr_df$t_1, dir_1)
  b <- pick(rr_df$t_2, dir_2)
  obs <- length(intersect(a, b))
  expected <- k^2 / n_shared
  # Each corner is tested in the direction it actually departs. Reporting only
  # the upper tail would print p = 1 on the discordant corners, hiding a
  # depletion significant at 1e-8 behind a number that reads as no result.
  enriched <- obs > expected
  list(
    k = k, observed = obs, expected = expected, enriched = enriched,
    direction = if (enriched) "enriched" else "depleted",
    p = if (enriched) {
      phyper(obs - 1, k, n_shared - k, k, lower.tail = FALSE)
    } else {
      phyper(obs, k, n_shared - k, k, lower.tail = TRUE)
    }
  )
}

quad_hyper <- list(
  UU = quadrant_hyper(1, 1), DD = quadrant_hyper(-1, -1),
  UD = quadrant_hyper(1, -1), DU = quadrant_hyper(-1, 1)
)

# The corner labels print a fold-change with no depth attached, so the depth is
# named here. The figure carries the rule and the caption carries the fact that
# fold-enrichment falls as the depth grows; three depths on the plot said the
# same thing in more ink.
sweep_txt <- sprintf("corners: top %.0f%% of each ranking", 100 * TOP_FRAC)

pE_subtitle <- if (descriptive) {
  sprintf(cfg$subtitle_fmt, n_shared)
} else {
  sprintf(cfg$subtitle_fmt, sweep_txt, n_shared)
}

message(sprintf(
  "  top-%d hypergeometric: %s",
  TOP_K,
  paste(sprintf(
    "%s obs=%d exp=%.1f %s p=%.3g", names(quad_hyper),
    vapply(quad_hyper, \(q) q$observed, numeric(1)),
    vapply(quad_hyper, \(q) q$expected, numeric(1)),
    vapply(quad_hyper, \(q) q$direction, character(1)),
    vapply(quad_hyper, \(q) q$p, numeric(1))
  ), collapse = " | ")
))

# stratified hypergeometric test (Cahill et al. 2018)
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_1, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_2, stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels = cfg$rrho_labels,
  log10.ind = TRUE,
  multipleTesting = "none",
  boundary = 0.02,
  method = "hyper",
  stepsize = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat)
nc <- ncol(hmat)
message(sprintf("  RRHO2 matrix: %d x %d", nr, nc))

na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

if (length(na_rows) && length(na_cols)) {
  row_before <- 1:(min(na_rows) - 1)
  row_after <- (max(na_rows) + 1):nr
  col_before <- 1:(min(na_cols) - 1)
  col_after <- (max(na_cols) + 1):nc
} else {
  mid <- floor(nr / 2)
  row_before <- 1:mid
  row_after <- (mid + 1):nr
  col_before <- 1:mid
  col_after <- (mid + 1):nc
}

# Both lists sorted descending (Up first):
# [row_before, col_before] = top of both lists = UU
# [row_after,  col_after]  = bottom of both    = DD
# [row_before, col_after]  = top list1, bottom list2 = UD
# [row_after,  col_before] = bottom list1, top list2 = DU
max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
max_DD <- max(hmat[row_after, col_after], na.rm = TRUE)
max_UD <- max(hmat[row_before, col_after], na.rm = TRUE)
max_DU <- max(hmat[row_after, col_before], na.rm = TRUE)

hotspot_genes <- list(
  UU = rrho_obj$genelist_uu$gene_list_overlap_uu,
  DD = rrho_obj$genelist_dd$gene_list_overlap_dd,
  UD = rrho_obj$genelist_ud$gene_list_overlap_ud,
  DU = rrho_obj$genelist_du$gene_list_overlap_du
)

n_UU <- length(hotspot_genes$UU)
n_DD <- length(hotspot_genes$DD)
n_UD <- length(hotspot_genes$UD)
n_DU <- length(hotspot_genes$DU)
message(sprintf("  Hotspot genes: UU=%d, DD=%d, UD=%d, DU=%d", n_UU, n_DD, n_UD, n_DU))

# no multiplier — prevents corner overlap in composite
txt_quad <- scale_text(BASE_QUADRANT, 146)

# jet colormap (Cahill et al. 2018 canonical appearance)
JET_COLORS <- c(
  "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
  "yellow", "#FF7F00", "red", "#7F0000"
)

# RRHO2 indexes from the up end of both lists, which puts up/up at the low-index
# corner and renders the panel as a mirror image of every published RRHO. Both
# axes are reflected so up/up sits top right, matching the convention readers
# arrive with (Plaisier et al. 2010; Cahill et al. 2018; Piron et al. 2023).
hmat_df <- expand.grid(row = 1:nr, col = 1:nc) |>
  mutate(
    neg_log10_p = as.vector(hmat),
    row = nr + 1L - row,
    col = nc + 1L - col
  )

max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)

# Anchors named by what the corner means rather than where it sits, so the
# reflection above cannot silently swap a label onto the wrong quadrant.
ann_x_up <- nr + 1L - min(row_before)
ann_x_down <- nr + 1L - max(row_after)
ann_y_up <- nc + 1L - min(col_before)
ann_y_down <- nc + 1L - max(col_after)

LABEL_FILL <- scales::alpha("white", 0.85)
LABEL_PADDING <- unit(0.45 * PRINT_SCALE, "mm")

# Each corner reports its hotspot size and the fixed-depth hypergeometric test,
# so the two empty corners read as a measured absence of overlap rather than as
# blank space. The heatmap summit is deliberately not quoted: it is a maximum
# over the rank grid, not a test statistic.
# Fold change rather than the raw counts. Two corner labels share each
# horizontal edge and grow toward each other, so "93 vs 30 enriched" collides
# where "3.1x enriched" does not. The counts stay in rrho2_summary.csv.
#
# Two rows, not four: a four-line box in each corner covered the heatmap it was
# annotating. The count rides on the title line because it names the quadrant as
# much as the words do; fold and test share the second line because they are one
# reading. The spaces come out of the p-value because the two labels sharing an
# edge have only about 280 pt between them, and at the spaced form their padded
# boxes touch.
quad_annotated <- function(text, n, hyper) {
  if (descriptive) {
    return(sprintf("%s\nn = %d", text, n))
  }
  fold <- if (hyper$enriched) {
    hyper$observed / hyper$expected
  } else {
    hyper$expected / max(hyper$observed, 0.5)
  }
  sprintf(
    "%s  n = %d\n%.1f× %s, %s",
    text, n, fold, hyper$direction, gsub(" ", "", fmt_p(hyper$p))
  )
}

ql <- list(
  UU = quad_annotated(cfg$quadrant_labels$UU, n_UU, quad_hyper$UU),
  DD = quad_annotated(cfg$quadrant_labels$DD, n_DD, quad_hyper$DD),
  UD = quad_annotated(cfg$quadrant_labels$UD, n_UD, quad_hyper$UD),
  DU = quad_annotated(cfg$quadrant_labels$DU, n_DU, quad_hyper$DU)
)

# geom_tile rather than geom_raster: a raster is one image, and the quartz PDF
# device draws rasterGrob vertically mirrored, which silently undid the row
# reflection above and left the four quadrant labels describing the wrong
# corners. Drawn as tiles the panel is vector, so no device can reorient it.
# The grid is stepsize 20, about 107 x 107 cells, so the cost is negligible.
pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = JET_COLORS,
    limits = c(0, max_val),
    na.value = "white",
    name = expression(-log[10](P)),
    guide = guide_colorbar(
      barwidth = unit(15 * PRINT_SCALE, "mm"), barheight = unit(2 * PRINT_SCALE, "mm"),
      title.position = "left", title.vjust = 0.5,
      title.theme = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey15")
    )
  ) +
  annotate("label",
    x = ann_x_up, y = ann_y_up,
    label = ql$UU,
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
    fontface = "bold", size = txt_quad,
    hjust = 1, vjust = 1
  ) +
  annotate("label",
    x = ann_x_down, y = ann_y_down,
    label = ql$DD,
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
    fontface = "bold", size = txt_quad,
    hjust = 0, vjust = 0
  ) +
  annotate("label",
    x = ann_x_up, y = ann_y_down,
    label = ql$UD,
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
    fontface = "bold", size = txt_quad,
    hjust = 1, vjust = 0
  ) +
  annotate("label",
    x = ann_x_down, y = ann_y_up,
    label = ql$DU,
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, label.r = unit(0.5, "mm"),
    fontface = "bold", size = txt_quad,
    hjust = 0, vjust = 1
  ) +
  scale_x_continuous(expand = expansion(mult = 0.015)) +
  scale_y_continuous(expand = expansion(mult = 0.015)) +
  labs(
    title = cfg$title,
    subtitle = pE_subtitle,
    x = cfg$axis_label_1,
    y = cfg$axis_label_2
  ) +
  FIG_THEME +
  theme(
    axis.text = element_blank(),
    axis.title.x = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 2)),
    axis.title.y = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(r = 2)),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = FIG_LEGEND_TEXT, face = "bold"),
    # coord_fixed makes the map as big as the smaller of its two free
    # dimensions, and that is the height: every point of dead space between the
    # axis title and the colour bar comes straight off the square.
    legend.margin = margin(0, 24, 0, 0, "mm"),
    legend.box.spacing = unit(1, "mm"),
    plot.margin = margin(0, 0, 0, 0, "mm")
  ) +
  coord_fixed(ratio = 1, clip = "off")

he <- cfg$hotspot_export_names
hotspot_export <- bind_rows(
  tibble(quadrant = he$UU, gene = hotspot_genes$UU),
  tibble(quadrant = he$DD, gene = hotspot_genes$DD),
  tibble(quadrant = he$UD, gene = hotspot_genes$UD),
  tibble(quadrant = he$DU, gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

ora_min <- cfg$ora_min_size %||% 15
pw_collection_E <- build_pathway_collection(
  min_size = ora_min, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)
all_genes_E <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (descriptive || length(gene_set) < 5) {
    return(tibble())
  }
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes_E,
      pathways       = pw_collection_E,
      em_cutoff      = 0.5,
      min_size       = ora_min,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) {
      message("  ORA error: ", e$message)
      tibble()
    }
  )
  if (nrow(res) > 0) {
    res |>
      mutate(
        quadrant = quadrant_name,
        pathway_label = clean_pathway_name(pathway),
        ID = pathway,
        p.adjust = padj,
        GeneRatio = paste0(overlap, "/", length(gene_set)),
        geneID = sapply(overlapGenes, paste, collapse = "/")
      ) |>
      arrange(padj, size)
  } else {
    tibble()
  }
}

oq <- cfg$ora_quadrant_names
ora_UU <- run_quadrant_ora(hotspot_genes$UU, oq$UU)
ora_DD <- run_quadrant_ora(hotspot_genes$DD, oq$DD)
ora_UD <- run_quadrant_ora(hotspot_genes$UD, oq$UD)
ora_DU <- run_quadrant_ora(hotspot_genes$DU, oq$DU)

og <- cfg$ora_grouped
ora_group_1 <- bind_rows(mget(og$file_1_quads))
ora_group_2 <- bind_rows(mget(og$file_2_quads))

# The note goes in the file itself, not beside it. As a sidecar it never
# reached the workbook, so the sheet shipped with no header and no rows and a
# reader could not tell "no enrichment" from "this step failed".
if (!is.null(og$note_if_empty_2) && nrow(ora_group_2) == 0) {
  ora_group_2 <- tibble(note = og$note_if_empty_2)
}

write_csv(ora_group_1, file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_group_2, file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

ggsave(file.path(RPT_PNG, paste0(file_stem, ".png")), pE_heat,
  width = PE_W, height = PE_W, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, paste0(file_stem, ".pdf")), pE_heat,
  width = PE_W, height = PE_W, units = "mm", device = pdf_device
)

MAX_PER_QUAD <- 12

if (!is.null(cfg$supp)) {
  txt_ora <- scale_text(BASE_STAT, PE_W)
  ora_all <- bind_rows(ora_group_1, ora_group_2)

  if (nrow(ora_all) > 0) {
    bar_df <- ora_all |>
      mutate(
        neg_log10_padj = -log10(p.adjust),
        pathway_label  = str_trunc(clean_pathway_name(pathway), 40),
        quadrant       = factor(quadrant, levels = cfg$supp$ora_quad_order)
      ) |>
      filter(!is.na(quadrant)) |>
      group_by(quadrant, .drop = FALSE) |>
      arrange(desc(neg_log10_padj)) |>
      slice_head(n = MAX_PER_QUAD) |>
      ungroup() |>
      filter(!is.na(neg_log10_padj)) |>
      arrange(quadrant, neg_log10_padj) |>
      mutate(uid = fct_inorder(paste0(pathway_label, "___", quadrant)))

    n_shown <- nrow(bar_df)
    n_total <- nrow(ora_all)

    pE_ora <- ggplot(bar_df, aes(x = neg_log10_padj, y = uid, fill = quadrant)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = overlap),
        hjust = -0.3, size = txt_ora * 0.7,
        color = "grey30"
      ) +
      scale_y_discrete(labels = function(x) str_remove(x, "___.*$")) +
      scale_fill_manual(values = cfg$ora_colors, guide = "none") +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      facet_grid(quadrant ~ .,
        scales = "free_y", space = "free_y",
        labeller = labeller(quadrant = cfg$supp$ora_quad_short)
      ) +
      labs(
        title = cfg$supp$ora_bar_title,
        subtitle = if (n_shown < n_total) {
          sprintf("Top %d per quadrant (%d terms total)", MAX_PER_QUAD, n_total)
        } else {
          sprintf("%d terms total", n_total)
        },
        x = expression(-log[10](p[adj])),
        y = NULL
      ) +
      FIG_THEME +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 8, hjust = 0.5),
        strip.text.y = element_text(size = 7, face = "bold", angle = 0),
        strip.background = element_rect(fill = "grey95", color = NA),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.margin = margin(2, 4, 2, 2, "mm")
      )

    pE_combined <- (pE_heat | pE_ora) + plot_layout(widths = c(1, 1.3))

    ggsave(file.path(cfg$supp$rpt_png, "SUPP_panel_E_rrho2_ora.png"), pE_combined,
      width = 400, height = 220, units = "mm", dpi = 300
    )
    ggsave(file.path(cfg$supp$rpt_pdf, "SUPP_panel_E_rrho2_ora.pdf"), pE_combined,
      width = 400, height = 220, units = "mm", device = pdf_device
    )
  }
}

sq <- cfg$summary_quadrant_names
quad_meta <- list(
  list(name = oq$UU, slug = sq$UU_slug, data = ora_UU, n_hot = n_UU),
  list(name = oq$DD, slug = sq$DD_slug, data = ora_DD, n_hot = n_DD),
  list(name = oq$UD, slug = sq$UD_slug, data = ora_UD, n_hot = n_UD),
  list(name = oq$DU, slug = sq$DU_slug, data = ora_DU, n_hot = n_DU)
)

for (qm in quad_meta) {
  if (nrow(qm$data) == 0) next
  q_df <- qm$data |>
    mutate(
      neg_log10_padj = -log10(padj),
      pathway_label = clean_pathway_name(pathway)
    ) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = MAX_PER_QUAD)
  message(sprintf("  %s: %d ORA pathways (%d hotspot genes)", qm$slug, nrow(q_df), qm$n_hot))
}

rrho2_meta <- tibble(
  quadrant = c(sq$UU, sq$DD, sq$UD, sq$DU),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  top_k = TOP_K,
  top_k_observed = vapply(quad_hyper[c("UU", "DD", "UD", "DU")], \(q) q$observed, numeric(1)),
  top_k_expected = round(vapply(quad_hyper[c("UU", "DD", "UD", "DU")], \(q) q$expected, numeric(1)), 1),
  top_k_direction = vapply(quad_hyper[c("UU", "DD", "UD", "DU")], \(q) q$direction, character(1)),
  top_k_p = signif(vapply(quad_hyper[c("UU", "DD", "UD", "DU")], \(q) q$p, numeric(1)), 3),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

pE_legend_grob <- cowplot::get_plot_component(pE_heat, "guide-box-bottom", return_all = FALSE)
if (!is.null(pE_legend_grob)) {
  pE_legend_plot <- cowplot::ggdraw(pE_legend_grob)
  ggsave(file.path(RPT_PNG, paste0(sub("_rrho2$", "", file_stem), "_legend.png")), pE_legend_plot,
    width = 60, height = 14, units = "mm", dpi = 300
  )
}

pE_title <- cfg$title
pE_legend <- NULL
pE_heat <- pE_heat +
  labs(title = NULL, subtitle = NULL, tag = NULL) +
  coord_fixed(ratio = 1, clip = "off")

message(sprintf("%s Panel E done", cfg$fig_id))
