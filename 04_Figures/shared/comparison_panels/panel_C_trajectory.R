# Shared Panel C — how much of the young training response survives in old,
# across every reasonable definition of a young responder.
#
# Training(Young) and Training(Old) share no group mean, the only such pairing
# in this four-group design, so this comparison needs no null calibration.
#
# Selecting proteins on the young statistic and then reporting old/young is
# regression to the mean: a protein enters the top set partly because its young
# estimate caught an upward noise draw, and that draw does not repeat in old.
# The bias is severe -- naively the FDR set retains 26% of its young magnitude,
# cross-fitted it retains 43%. Every set is therefore cross-fitted: half the
# young subjects choose the members, the other half supply the young magnitude,
# so membership carries no information about the number being measured. The
# whole proteome needs no correction and anchors the panel.

source("04_Figures/shared/style.R")
source("04_Figures/shared/print_scale_apply.R")

pacman::p_load(dplyr, tibble, tidyr, ggplot2, forcats)

set.seed(42)

N_SPLIT <- 300L
FDR_CUT <- 0.05
BR_INSET <- 0.07

meta <- cfg$meta
meta$era <- sub("_(Pre|Post)$", "", meta$Group_Time)
E <- cfg$matrix

dep <- cfg$dep_df |> distinct(uniprot_id, .keep_all = TRUE)
E <- E[rownames(E) %in% dep$uniprot_id, , drop = FALSE]
dep <- dep[match(rownames(E), dep$uniprot_id), ]

# Subject_ID is not unique within an era: O_S06 and OP_S06 are two different
# people, as are O_S17 and OP_S17. Keying on it merged each pair into one
# subject carrying two Pre and two Post, which then failed the test below and
# took four old participants with usable biopsies out of the panel. The Col_ID
# stem is the identifier 00_input/README.md says to use.
subject_key <- function(col_id) sub("_(Pre|Post)$", "", col_id)

# Only subjects with exactly one Pre and one Post inside their own era are used,
# so no subject's Post is ever paired with another's Pre.
clean_subjects <- function(era) {
  m <- meta[meta$era == era, ]
  tb <- table(subject_key(m$Col_ID), sub(".*_", "", m$Group_Time))
  rownames(tb)[tb[, "Pre"] == 1 & tb[, "Post"] == 1]
}

paired_delta <- function(subj, era) {
  idx <- function(tp) {
    r <- meta[meta$Group_Time == paste0(era, "_", tp), ]
    r$Col_ID[match(subj, subject_key(r$Col_ID))]
  }
  rowMeans(E[, idx("Post"), drop = FALSE] - E[, idx("Pre"), drop = FALSE], na.rm = TRUE)
}

young_subj <- clean_subjects("Young")
old_subj <- clean_subjects("Old")
d_old <- paired_delta(old_subj, "Old")
d_young <- paired_delta(young_subj, "Young")

n_fdr <- sum(dep$adj.P.Val_Training_Young < FDR_CUT, na.rm = TRUE)
# sig_pi is -1/0/1. Testing it against 1 counted only the up-regulated half
# and sized this set at 43 while labelling it with the 99 the threshold gives.
n_pi <- sum(dep$sig_pi_Training_Young != 0, na.rm = TRUE)

# Each rule takes a young effect vector, so the same rule applies to the full
# cohort or to half of it without being restated.
SETS <- list(
  `All\nmeasured` = function(dy) rep(TRUE, length(dy)),
  `Top 500` = function(dy) rank(-abs(dy), ties.method = "first") <= 500,
  `Top 250` = function(dy) rank(-abs(dy), ties.method = "first") <= 250,
  `FDR < 0.05` = function(dy) rank(-abs(dy), ties.method = "first") <= n_fdr,
  `Π < 0.05` = function(dy) rank(-abs(dy), ties.method = "first") <= n_pi
)

splits <- replicate(N_SPLIT, {
  a <- sample(young_subj, floor(length(young_subj) / 2))
  b <- setdiff(young_subj, a)
  da <- paired_delta(a, "Young")
  db <- paired_delta(b, "Young")
  vapply(SETS, function(f) {
    k <- f(da)
    c(median(abs(db[k]), na.rm = TRUE), median(abs(d_old[k]), na.rm = TRUE))
  }, numeric(2))
})

naive_ratio <- vapply(SETS, function(f) {
  k <- f(d_young)
  median(abs(d_old[k]), na.rm = TRUE) / median(abs(d_young[k]), na.rm = TRUE)
}, numeric(1))

trajectory_stats <- bind_rows(lapply(seq_along(SETS), function(i) {
  y <- splits[1, i, ]
  o <- splits[2, i, ]
  tibble(
    set = names(SETS)[i],
    n = sum(SETS[[i]](d_young)),
    young_magnitude = median(y), old_magnitude = median(o),
    young_lo = quantile(y, 0.025), young_hi = quantile(y, 0.975),
    old_lo = quantile(o, 0.025), old_hi = quantile(o, 0.975),
    retained = median(o / y),
    retained_lo = quantile(o / y, 0.025),
    retained_hi = quantile(o / y, 0.975),
    retained_naive = naive_ratio[[i]]
  )
})) |>
  mutate(
    set = fct_inorder(set),
    # The bracket stops short of both bar tops by a fraction of what it spans,
    # so the wide groups read less heavily than the narrow ones. Proportional
    # rather than fixed: a fixed inset would swallow the All-measured gap,
    # which is a seventh of the Pi one.
    br_lo = old_magnitude + BR_INSET * (young_magnitude - old_magnitude),
    br_hi = young_magnitude - BR_INSET * (young_magnitude - old_magnitude)
  )

write.csv(
  trajectory_stats |> mutate(set = gsub("\n", " ", as.character(set))),
  file.path(cfg$dat, "panel_C_trajectory.csv"),
  row.names = FALSE
)

write.csv(
  tibble(
    gene = dep$gene,
    logFC_Training_Young = d_young,
    logFC_Training_Old = d_old,
    in_top_500 = SETS[["Top 500"]](d_young),
    in_top_250 = SETS[["Top 250"]](d_young),
    in_fdr = SETS[["FDR < 0.05"]](d_young),
    in_pi = SETS[["Π < 0.05"]](d_young)
  ) |> arrange(desc(abs(logFC_Training_Young))),
  file.path(cfg$dat, "panel_C_trajectory_proteins.csv"),
  row.names = FALSE
)

bars <- trajectory_stats |>
  select(set, young_magnitude, old_magnitude, young_lo, young_hi, old_lo, old_hi) |>
  pivot_longer(-set,
    names_to = c("era", ".value"),
    names_pattern = "^(young|old)_(.*)$"
  ) |>
  rename(mag = magnitude) |>
  mutate(era = factor(
    ifelse(era == "young", "Young", "Old"),
    levels = c("Young", "Old")
  ))

dodge <- position_dodge(width = 0.62)
BR_X <- 0.045
BR_CAP <- 0.045
BR_TICK <- 0.065
BR_LW <- 0.6
BR_COL <- "grey10"
y_top <- max(trajectory_stats$young_hi, trajectory_stats$old_hi)

pC_trajectory <- ggplot(bars, aes(set, mag, fill = era)) +
  geom_col(
    position = dodge, width = 0.58, colour = "black", linewidth = 0.2
  ) +
  geom_linerange(aes(ymin = lo, ymax = hi),
    position = dodge, linewidth = 0.35, colour = "grey25"
  ) +
  # A capped bracket rather than a diagonal, matching the abstract card. It sits
  # clear of the bars because the dodge leaves only 0.02 units between them,
  # about 0.4 mm, where a spine is invisible and the caps merge into the bar
  # edges. BR_X puts it just right of the old bar; the caps are short, so the
  # bracket reads by height alignment the way the abstract's does. Both caps
  # point right, making the shape a bracket rather than the capped error bars
  # already on every bar.
  geom_segment(
    data = trajectory_stats, inherit.aes = FALSE,
    aes(
      x = as.numeric(set) + BR_X, xend = as.numeric(set) + BR_X,
      y = br_lo, yend = br_hi
    ),
    colour = BR_COL, linewidth = BR_LW
  ) +
  geom_segment(
    data = trajectory_stats, inherit.aes = FALSE,
    aes(
      x = as.numeric(set) + BR_X - BR_CAP, xend = as.numeric(set) + BR_X,
      y = br_hi, yend = br_hi
    ),
    colour = BR_COL, linewidth = BR_LW
  ) +
  geom_segment(
    data = trajectory_stats, inherit.aes = FALSE,
    aes(
      x = as.numeric(set) + BR_X - BR_CAP, xend = as.numeric(set) + BR_X,
      y = br_lo, yend = br_lo
    ),
    colour = BR_COL, linewidth = BR_LW
  ) +
  # Each bar carries its own magnitude, centred inside it: above the bar the
  # label sits past the split-half interval and drifts away from what it names.
  geom_text(
    data = bars, inherit.aes = FALSE,
    aes(x = set, y = mag / 2, label = sub("^0", "", sprintf("%.2f", mag)), group = era),
    position = dodge, vjust = 0.5, size = (FIG_AXIS_TEXT - 4) / .pt,
    fontface = "bold", colour = "white"
  ) +
  # The retained fraction is what the panel is for, so it is set one point
  # above the axis type rather than off BASE_STAT, which print_scale_apply
  # deliberately holds back for panel A's fitted label boxes.
  geom_segment(
    data = trajectory_stats, inherit.aes = FALSE,
    aes(
      x = as.numeric(set) + BR_X, xend = as.numeric(set) + BR_X + BR_TICK,
      y = (young_magnitude + old_magnitude) / 2,
      yend = (young_magnitude + old_magnitude) / 2
    ),
    colour = BR_COL, linewidth = BR_LW
  ) +
  geom_text(
    data = trajectory_stats, inherit.aes = FALSE,
    aes(
      x = as.numeric(set) + BR_X + BR_TICK + 0.03,
      y = (young_magnitude + old_magnitude) / 2,
      # Stated as the fraction lost, not kept, to match the abstract card:
      # one number, one direction, in both places a reader meets it.
      label = sprintf("%.0f%%", 100 * (1 - retained))
    ),
    hjust = 0, vjust = 0.5, size = (FIG_AXIS_TEXT - 1) / .pt,
    fontface = "bold", colour = "grey15"
  ) +
  # The two series are the Training(Young) and Training(Old) contrasts, not age
  # groups, so they take the contrast palette F02 already uses for them. The
  # direction palette they used to borrow read blue for Young, which collides
  # with the Down colour everywhere else in the tree.
  scale_fill_manual(
    values = c(
      Young = CONTRAST_COLORS[["Training_Young"]],
      Old = CONTRAST_COLORS[["Training_Old"]]
    ),
    name = NULL
  ) +
  # Five categories share about 245 pt of axis, which one line of the
  # print-scaled type no longer fits. Wrapping here and not in the names keeps
  # the exports and the SETS lookups on the flat strings.
  # The percentage sits right of the bracket, so the last category needs room
  # past its own centre.
  scale_x_discrete(
    labels = \(x) sub(" ", "\n", x),
    expand = expansion(add = c(0.55, 0.70))
  ) +
  scale_y_continuous(expand = expansion(c(0, 0.20))) +
  labs(x = NULL, y = expression(bold("median |log"[2] * "FC|"))) +
  FIG_THEME +
  theme(
    axis.title.y = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    axis.text.x = element_text(
      size = FIG_AXIS_TEXT - 2, face = "bold", lineheight = 0.9
    ),
    # No rules behind the bars, as on the abstract card. Each bar prints its
    # own magnitude, so the horizontal gridlines were reading nothing.
    panel.grid = element_blank(),
    # Inside the plot box, where the tallest bars are not: the panel used to
    # reserve 52 pt below itself to offset a legend sitting above the plot, and
    # that reserve was what left C shorter than its row-mates.
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    # An in-plot position defaults to a vertical key stack, whose two labels
    # then sit close enough to touch.
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.size = unit(2.6, "mm"),
    legend.text = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    legend.margin = margin(0, 0, 0, 0),
    # Set so both borders land level with panel B's tile field in the F04
    # composite, measured with tools/measure_panels.py. C fills its cell, so
    # the bottom margin lifts the panel and its wrapped axis labels together;
    # matching B costs C about a tenth of its height, which is the trade the
    # alignment asks for.
    plot.margin = margin(8.4, 5.5, 35.6, 5.5)
  )

ggsave(file.path(cfg$rpt_png, "MAIN_panel_C_trajectory.png"), pC_trajectory,
  width = cfg$panel_w, height = cfg$panel_h, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(cfg$rpt_pdf, "MAIN_panel_C_trajectory.pdf"), pC_trajectory,
  width = cfg$panel_w, height = cfg$panel_h, units = "mm", device = get_pdf_device()
)

message(sprintf(
  "Panel C: %d young / %d old paired subjects | %s",
  length(young_subj), length(old_subj),
  paste(sprintf(
    "%s %.0f%% (naive %.0f%%)",
    gsub("\n", " ", as.character(trajectory_stats$set)),
    100 * trajectory_stats$retained, 100 * trajectory_stats$retained_naive
  ), collapse = " | ")
))
