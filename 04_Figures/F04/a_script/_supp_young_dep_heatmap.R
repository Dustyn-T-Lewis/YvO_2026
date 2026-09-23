#!/usr/bin/env Rscript
# F04 Supplementary (S9 Figure): protein-level training response across age
# groups. The rest of the paper reports the training response as model
# coefficients and pathway scores; this is the only view of per-protein
# abundance across individual participants. Rows are the proteins reaching
# FDR < 0.05 with training in younger adults, split by whether the older-adult
# response shares the sign of the younger one, with over-representation
# analysis of the concordant rows.
#
# Display is the paired delta (Post - Pre) per subject, so the two younger
# participants without both timepoints (Y_S07 Pre only, Y_S05 Post only) drop
# out and 15 Young + 15 Old subjects remain. The limma logFC annotation beside
# the matrix is estimated on all 62 samples and is unaffected by that.
#
# Sourced by 02_supp_panels.R — expects style.R already loaded.

source("04_Figures/shared/pathway_utils.R")

pacman::p_load(ComplexHeatmap, circlize, gridExtra)

BASE <- "04_Figures/F04"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT <- file.path(BASE, "c_data", "panel_supp")
for (d in c(RPT_PNG, RPT_PDF, DAT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

pdf_device <- get_pdf_device()

dep <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

id_cols <- c("uniprot_id", "protein", "gene", "description")
sample_cols <- setdiff(names(dep), id_cols)
sample_cols <- sample_cols[grepl("_(Pre|Post)$", sample_cols)]
subject_key <- unique(sub("_(Pre|Post)$", "", sample_cols))
paired_subj <- subject_key[
  paste0(subject_key, "_Pre") %in% sample_cols &
    paste0(subject_key, "_Post") %in% sample_cols
]
unpaired <- setdiff(subject_key, paired_subj)
paired_subj <- paired_subj[order(!grepl("^Y", paired_subj), paired_subj)]

young_dep <- dep |>
  filter(!is.na(adj.P.Val_Training_Young), adj.P.Val_Training_Young < 0.05) |>
  mutate(response = if_else(
    sign(logFC_Training_Young) == sign(logFC_Training_Old),
    "Concordant", "Discordant"
  )) |>
  mutate(direction = if_else(logFC_Training_Young > 0, "Up", "Down")) |>
  arrange(factor(direction, levels = c("Up", "Down")), desc(logFC_Training_Young))

delta <- vapply(
  paired_subj,
  function(s) young_dep[[paste0(s, "_Post")]] - young_dep[[paste0(s, "_Pre")]],
  numeric(nrow(young_dep))
)
rownames(delta) <- young_dep$gene

age_of_subj <- if_else(grepl("^Y", paired_subj), "Young", "Old")
col_split <- factor(age_of_subj, levels = c("Young", "Old"))
row_split <- factor(young_dep$direction, levels = c("Up", "Down"))

cap <- unname(quantile(abs(delta), 0.98, na.rm = TRUE))
col_fun <- colorRamp2(c(-cap, 0, cap), c("#4393C3", "white", "#D6604D"))
lfc_cap <- max(abs(c(
  young_dep$logFC_Training_Young,
  young_dep$logFC_Training_Old
)))
lfc_fun <- colorRamp2(c(-lfc_cap, 0, lfc_cap), c("#4393C3", "white", "#D6604D"))

CONC_COLORS <- c(Concordant = "#5AAE61", Discordant = "#9970AB")

ha_row <- rowAnnotation(
  Concordance = young_dep$response,
  `logFC Young` = young_dep$logFC_Training_Young,
  `logFC Old` = young_dep$logFC_Training_Old,
  col = list(
    Concordance = CONC_COLORS,
    `logFC Young` = lfc_fun, `logFC Old` = lfc_fun
  ),
  annotation_name_gp = gpar(fontsize = 7),
  annotation_name_rot = 45,
  show_legend = c(Concordance = TRUE, `logFC Young` = TRUE, `logFC Old` = FALSE),
  annotation_legend_param = list(
    Concordance = list(
      title = "Older-adult response",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    ),
    `logFC Young` = list(
      title = "limma logFC", direction = "horizontal",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    )
  )
)

ht <- Heatmap(
  delta,
  name = "delta_log2", col = col_fun, na_col = "grey92",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_split = row_split, column_split = col_split,
  row_gap = unit(2, "mm"), column_gap = unit(3, "mm"),
  # 135 rows at column width leaves 2.17 pt per gene name, a third the height
  # of body text and unreadable in print. The identities are in S5 Table and
  # the same protein set is drawn with legible labels in S4b, so the rows are
  # unlabelled here and the panel carries the pattern rather than the names.
  show_row_names = FALSE,
  column_names_gp = gpar(fontsize = 6),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  right_annotation = ha_row,
  heatmap_legend_param = list(
    title = "log2 delta (Post - Pre)",
    direction = "horizontal",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  )
)

n_conc <- sum(young_dep$response == "Concordant")
n_disc <- sum(young_dep$response == "Discordant")
rho <- cor(young_dep$logFC_Training_Young, young_dep$logFC_Training_Old,
  method = "spearman"
)

fig_w_mm <- 320
fig_h_mm <- max(180, 1.9 * nrow(delta) + 60)
title_grob <- textGrob(
  sprintf(
    "Protein-level training response across age groups  |  %d younger-adult FDR DEPs  |  %d concordant, %d discordant  |  Spearman rho = %.3f",
    nrow(delta), n_conc, n_disc, rho
  ),
  gp = gpar(fontsize = 8, fontface = "bold")
)

g_ht <- grid.grabExpr(draw(ht,
  heatmap_legend_side = "bottom",
  merge_legend = TRUE,
  padding = unit(c(2, 4, 2, 14), "mm")
))

pw_list <- build_pathway_collection(
  min_size = 10, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)
universe <- unique(dep$gene[!is.na(dep$gene)])

# The flow panel needs pw_list and universe, and contributes the right-hand
# third of the composite.
source("04_Figures/F04/a_script/_supp_young_dep_ora_flow.R")
g_bars <- ggplotGrob(pS_ora_bars)

write_composite <- function() {
  grid.arrange(g_ht, g_bars, ncol = 2, widths = c(2.1, 1), top = title_grob)
}

png(file.path(RPT_PNG, "SUPP_young_dep_heatmap.png"),
  width = fig_w_mm, height = fig_h_mm, units = "mm", res = 300
)
write_composite()
dev.off()

pdf_dev <- get_pdf_device()
if (is.character(pdf_dev)) pdf_dev <- match.fun(pdf_dev)
pdf_dev(file.path(RPT_PDF, "SUPP_young_dep_heatmap.pdf"),
  width = fig_w_mm / 25.4, height = fig_h_mm / 25.4
)
write_composite()
dev.off()

set.seed(42)
ora_conc <- run_ora_deduplicated(
  genes = young_dep$gene[young_dep$response == "Concordant"],
  universe = universe, pathways = pw_list,
  em_cutoff = 0.5, min_size = 10, max_size = 500
)

# The 25 discordant proteins are too few for over-representation testing, so
# they are exported per protein rather than summarised as pathways.
disc_df <- young_dep |>
  filter(response == "Discordant") |>
  transmute(gene, description,
    logFC_Training_Young = round(logFC_Training_Young, 3),
    logFC_Training_Old = round(logFC_Training_Old, 3),
    adj.P.Val_Training_Young = signif(adj.P.Val_Training_Young, 3),
    adj.P.Val_Training_Old = signif(adj.P.Val_Training_Old, 3)
  ) |>
  arrange(desc(abs(logFC_Training_Young)))

heat_df <- as.data.frame(delta) |>
  mutate(
    gene = rownames(delta), response = as.character(row_split),
    .before = 1
  )

write_csv(heat_df, file.path(DAT, "SUPP_young_dep_heatmap.csv"))
write_csv(ora_conc, file.path(DAT, "SUPP_young_dep_ora_concordant.csv"))
write_csv(disc_df, file.path(DAT, "SUPP_young_dep_discordant.csv"))

message(sprintf(
  "F04 young-DEP heatmap: %d proteins x %d paired subjects (%s unpaired dropped), %d concordant, %d discordant, %d ORA terms",
  nrow(delta), ncol(delta), paste(unpaired, collapse = "/"),
  n_conc, n_disc, nrow(ora_conc)
))
