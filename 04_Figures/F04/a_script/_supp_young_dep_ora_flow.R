#!/usr/bin/env Rscript
# F04 Supplementary: what the younger-adult FDR set is, beside the heatmap that
# shows how it moves.
#
# Over-representation runs separately on the increases and the decreases. The
# two enrich for opposite biology, so testing them together dilutes both: one
# test on all 135 returns 28 terms, split it returns 52 and 13.
#
# Expects: young_dep (gene, logFC_Training_Young, direction) plus universe,
# pw_list, DAT and style.R. Exports pS_ora_bars.

N_TERMS <- 5L
DIR_FILL <- c(Increased = "#C0483B", Decreased = "#2E6F9E")

short_term <- function(x) {
  clean_pathway_name(x) |>
    stringr::str_remove("^Reference ") |>
    stringr::str_replace("Microtubule", "MT") |>
    stringr::str_replace("Electron Transport", "ETC") |>
    stringr::str_replace("Mitochondrial", "Mito.") |>
    stringr::str_replace("Ubiquinone", "UQ") |>
    stringr::str_replace("Transmembrane", "Transmem.") |>
    stringr::str_replace("Copi Independent", "Coat Independent") |>
    stringr::str_replace("Retrograde Traffic", "Retro. Trans.") |>
    stringr::str_replace("Intermediates By Cct Tric", "Intermediates") |>
    stringr::str_replace(" And ", " & ") |>
    stringr::str_wrap(width = 30)
}

ora_direction <- function(dir_label) {
  genes <- young_dep$gene[young_dep$direction == dir_label]
  set.seed(42)
  run_ora_deduplicated(
    genes = genes, universe = universe, pathways = pw_list,
    em_cutoff = 0.5, min_size = 10, max_size = 500
  ) |>
    filter(dedup_status == "kept", padj < 0.05) |>
    arrange(padj) |>
    head(N_TERMS) |>
    mutate(group = if (dir_label == "Up") "Increased" else "Decreased", n_set = length(genes))
}

ora_tbl <- bind_rows(lapply(c("Up", "Down"), ora_direction)) |>
  mutate(group = factor(group, levels = c("Increased", "Decreased")))

write_csv(
  ora_tbl |>
    mutate(term = clean_pathway_name(pathway), neg_log10_fdr = -log10(padj)) |>
    select(group, term, database_term = pathway, overlap, size, padj, neg_log10_fdr),
  file.path(DAT, "SUPP_young_dep_ora.csv")
)

bar_df <- ora_tbl |>
  mutate(neg_log10 = -log10(padj)) |>
  arrange(group, desc(neg_log10)) |>
  mutate(
    slot = row_number() + (as.integer(group) - 1) * 2,
    y = max(slot) - slot + 1,
    lab = sprintf("%s  (%d/%d)", short_term(pathway), overlap, size)
  )

head_df <- bar_df |>
  group_by(group, n_set) |>
  summarise(y = max(y) + 1.1, .groups = "drop") |>
  mutate(lab = sprintf("%s with training in younger adults (n = %d)", group, n_set))

x_max <- max(bar_df$neg_log10)

pS_ora_bars <- ggplot(bar_df, aes(x = neg_log10, y = y, fill = group)) +
  geom_col(width = 0.62, orientation = "y") +
  geom_text(aes(x = -0.3, label = lab),
    hjust = 1, size = 2.1, lineheight = 0.95, colour = "grey15"
  ) +
  geom_text(
    data = head_df, aes(x = -0.3, y = y, label = lab),
    inherit.aes = FALSE, hjust = 1, size = 2.3, fontface = "bold", colour = "grey15"
  ) +
  scale_fill_manual(values = DIR_FILL, guide = "none") +
  scale_x_continuous(
    limits = c(-x_max * 1.45, x_max * 1.06),
    breaks = scales::pretty_breaks(4)(c(0, x_max)), expand = c(0, 0)
  ) +
  scale_y_continuous(expand = expansion(add = 1.2)) +
  labs(title = "Enriched pathways", x = expression(-log[10] ~ FDR), y = NULL) +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 7),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 6.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.2, colour = "grey88"),
    plot.margin = margin(2, 3, 2, 1, "mm")
  )

message(sprintf("  ORA: %d terms over %d directions", nrow(bar_df), 2L))
