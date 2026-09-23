# Shared Panel F — Training(Young) signature positions in the Training(Old)
# ranking, as a barcode with enrichment worms.
#
# This is limma's barcodeplot idiom (Wu & Smyth 2012, NAR 40:e133) drawn in
# ggplot so it composes into the figure. Up-regulated members sit above the
# axis, down-regulated below, and the worm is a tricube moving average of set
# membership along the rank, so a rise means members are concentrated there.
#
# Only Training(Young)-defined sets appear. They are the one selection this
# design can test against the Training(Old) ranking without calibration:
# Young_Post - Young_Pre and Old_Post - Old_Pre share no group mean. Sets
# defined by Aging share Old_Pre with the ranking and need a swap null instead,
# which is a different argument and belongs elsewhere.

source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_apply.R")

pacman::p_load(dplyr, tibble, tidyr, ggplot2, limma)

set.seed(42)

WORM_SPAN <- 0.12
TICK_INNER <- 0.06
TICK_OUTER <- 0.30
WORM_FLOOR <- 0.34

dep <- cfg$dep_df
abundance <- cfg$matrix
subject <- sub("_(Pre|Post)$", "", cfg$meta$Col_ID)

design <- model.matrix(~ 0 + cfg$group_time)
colnames(design) <- levels(cfg$group_time)
corfit <- duplicateCorrelation(abundance, design, block = subject)
tro <- makeContrasts(Training_Old = Old_Post - Old_Pre, levels = design)

ranking <- dep$t_Training_Old
names(ranking) <- dep$uniprot_id
ranking <- ranking[rownames(abundance)]
order_by_rank <- order(ranking, decreasing = TRUE)
n_ranked <- length(order_by_rank)

member_rows <- function(criterion, direction) {
  col <- paste0(criterion, "_Training_Young")
  lfc <- dep$logFC_Training_Young
  hit <- !is.na(dep[[col]]) & dep[[col]] < 0.05 &
    if (direction == "up") lfc > 0 else lfc < 0
  which(rownames(abundance) %in% dep$uniprot_id[hit])
}

CRITERIA <- c(FDR = "adj.P.Val", `Π` = "pi_score")
sets <- list()
for (crit in names(CRITERIA)) {
  for (dir in c("up", "down")) {
    sets[[paste(crit, dir)]] <- member_rows(CRITERIA[[crit]], dir)
  }
}

# What each set is made of. fry says whether a set sits high in the older
# ranking; it says nothing about what the set does, and the same four sets carry
# that answer already. Collection settings match the retired panel_C_fry.R, so
# these terms are comparable with the ones that figure reported.
source("04_Figures/shared/pathway_utils.R")

ora_pw <- build_pathway_collection(
  min_size = 10, max_size = 500, include_goslim = TRUE, exclude_variants = TRUE
)

top_term <- function(rows) {
  genes <- unique(dep$gene[match(rownames(abundance)[rows], dep$uniprot_id)])
  genes <- genes[!is.na(genes)]
  if (length(genes) < 5) {
    return(NA_character_)
  }
  r <- run_ora_deduplicated(
    genes = genes, universe = unique(dep$gene[!is.na(dep$gene)]),
    pathways = ora_pw, min_size = 10, max_size = 500, padj_cutoff = 0.1
  )
  if (!nrow(r)) {
    return(NA_character_)
  }
  # GO Slim and Hallmark first. Both are curated summary layers, so their names
  # say what a reader needs in a few words, where the Reactome winner for the
  # same genes was "Rab7 Regulated Microtubule Minus End Directed Transport".
  # Everything else only if neither has a term for this set.
  summary_hits <- r[r$database %in% c("GO Slim", "Hallmark"), ]
  if (nrow(summary_hits)) {
    r <- summary_hits
  }
  # "Reference " prefixes a KEGG family name and carries nothing, and the
  # articles carry little. A lookahead, not a plain match: " Of The " is two
  # overlapping matches sharing one space, and a plain gsub leaves the second.
  # Title case mangles the RNA species, which the collection writes as one word.
  tidy <- function(x) {
    clean_pathway_name(x) |>
      sub(pattern = "^Reference ", replacement = "") |>
      gsub(pattern = " (Of|The|And)(?= )", replacement = "", perl = TRUE) |>
      gsub(pattern = "\\b([TMR])rna\\b", replacement = "\\L\\1\\ERNA", perl = TRUE)
  }
  # Terms within a factor of two of the best are not meaningfully separable at
  # this sample size, so the shortest name wins rather than the longest.
  near <- r[r$padj <= 2 * min(r$padj), ]
  near$label <- vapply(near$pathway, tidy, character(1))
  r <- near[which.min(nchar(near$label)), ]
  term <- r$label
  # The six-column box runs out near 34 characters; cut at a word, not through.
  if (nchar(term) > 34) {
    term <- sub(" [^ ]*$", "", strtrim(term, 35))
  }
  sprintf("%s, %s", term, gsub(" ", "", fmt_p(r$padj)))
}

ora_top <- vapply(sets, top_term, character(1))

fry_res <- fry(abundance,
  index = sets[lengths(sets) >= 3], design = design,
  contrast = tro[, "Training_Old"],
  block = subject, correlation = corfit$consensus.correlation
)
fry_p <- setNames(fry_res$PValue, rownames(fry_res))

# Membership along the rank, smoothed the way limma's barcodeplot does it.
worm_of <- function(rows) {
  membership <- as.integer(seq_len(n_ranked) %in% match(rows, order_by_rank))
  tricubeMovingAverage(membership, span = WORM_SPAN)
}

lane <- function(crit, dir) {
  key <- paste(crit, dir)
  rows <- sets[[key]]
  positions <- match(rows, order_by_rank)
  worm <- worm_of(rows)
  side <- if (dir == "up") 1 else -1
  scaled <- WORM_FLOOR + (1 - WORM_FLOOR) * worm / max(worm)
  list(
    ticks = tibble(
      criterion = crit, rank = positions,
      y0 = side * TICK_INNER, y1 = side * TICK_OUTER
    ),
    worm = tibble(
      criterion = crit, rank = seq_len(n_ranked),
      y = side * scaled, y0 = side * WORM_FLOOR, direction = dir
    ),
    # Each worm peaks at its own end of the ranking, so the label goes to the
    # opposite end and stays off the curve.
    label = tibble(
      criterion = crit, direction = dir,
      x = if (dir == "up") n_ranked * 0.995 else n_ranked * 0.005,
      hjust = if (dir == "up") 1 else 0,
      # Anchored to the corner, not centred on it: vjust lets a two-line block
      # grow down from the top on the up lane and up from the floor on the
      # down lane, so neither runs off its own track.
      y = if (dir == "up") 0.97 else -0.97,
      vjust = if (dir == "up") 1 else 0,
      text = paste(
        sprintf("%s (n = %d), %s", dir, length(rows), fmt_p(fry_p[[key]])),
        ora_top[[key]],
        sep = "\n"
      )
    )
  )
}

lanes <- lapply(names(CRITERIA), \(crit) lapply(c("up", "down"), \(d) lane(crit, d)))
lanes <- unlist(lanes, recursive = FALSE)
pull_part <- function(part) bind_rows(lapply(lanes, `[[`, part))

ticks <- pull_part("ticks")
worms <- pull_part("worm")
labels <- pull_part("label")

criterion_levels <- names(CRITERIA)
as_facet <- function(d) mutate(d, criterion = factor(criterion, levels = criterion_levels))

fry_barcode_stats <- list(
  n_ranked = n_ranked,
  sets = lengths(sets),
  p_values = fry_p[names(sets)]
)

write.csv(
  tibble(set = names(sets), n = lengths(sets), fry_p = fry_p[names(sets)]),
  file.path(cfg$dat, "panel_F_fry_barcode.csv"),
  row.names = FALSE
)

pF_barcode <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.35) +
  geom_segment(
    data = as_facet(ticks),
    aes(x = rank, xend = rank, y = y0, yend = y1, colour = y1 > 0),
    linewidth = 0.2, alpha = 0.8
  ) +
  geom_ribbon(
    data = as_facet(worms),
    aes(rank, ymin = y0, ymax = y, fill = direction == "up"), alpha = 0.25
  ) +
  geom_line(
    data = as_facet(worms),
    aes(rank, y, colour = direction == "up"), linewidth = 0.5
  ) +
  geom_text(
    data = as_facet(labels),
    aes(x, y, label = text, hjust = hjust, vjust = vjust, colour = y > 0),
    # These four strings are the panel's result, so they match the axis type
    # rather than BASE_STAT, which print_scale_apply holds at the older factor.
    size = FIG_AXIS_TEXT / .pt, fontface = "bold.italic", show.legend = FALSE
  ) +
  facet_wrap(~criterion, ncol = 1, strip.position = "left") +
  scale_colour_manual(
    values = c(`TRUE` = DIR_COLORS[["Up"]], `FALSE` = DIR_COLORS[["Down"]]),
    aesthetics = c("colour", "fill"), guide = "none"
  ) +
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(limits = c(-1.05, 1.05), breaks = NULL, expand = expansion(0)) +
  labs(
    x = sprintf("Training (Old) rank, up to down  (n = %d)", n_ranked),
    y = "Training (Young) set members"
  ) +
  FIG_THEME +
  theme(
    strip.placement = "outside", strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 90, size = FIG_AXIS_TEXT + 1, face = "bold"),
    axis.title.x = element_text(
      size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 2)
    ),
    axis.title.y = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    axis.text.x = element_text(size = FIG_AXIS_TEXT - 0.5),
    panel.grid = element_blank(), panel.spacing.y = unit(1.5, "mm"),
    # Both margins trimmed to sit level with the NES scatters, measured against
    # their border in the rendered composite. The axis furniture below the
    # facets is a fixed 32 pt, so dropping the top margin raises the top edge by
    # exactly 5.5 and leaves the bottom where it is.
    plot.margin = margin(0, 5.5, 0, 5.5)
  )

ggsave(file.path(cfg$rpt_png, "MAIN_panel_F_fry_barcode.png"), pF_barcode,
  width = cfg$panel_w, height = cfg$panel_h, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(cfg$rpt_pdf, "MAIN_panel_F_fry_barcode.pdf"), pF_barcode,
  width = cfg$panel_w, height = cfg$panel_h, units = "mm", device = get_pdf_device()
)

message(sprintf(
  "Panel F: %d ranked proteins | %s",
  n_ranked,
  paste(sprintf("%s n=%d %s", names(sets), lengths(sets), vapply(fry_p[names(sets)], fmt_p, character(1))),
    collapse = " | "
  )
))
