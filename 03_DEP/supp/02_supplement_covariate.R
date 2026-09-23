#!/usr/bin/env Rscript
# R1.2: "Did the authors consider including 'Parent Study Supplementation
# Group' as a random effect or covariate in their empirical Bayes limma model?"
#
# Four questions, answered against the design matrix rather than by argument.
#
#   1. Can supplement enter the model that estimates Aging? Not usefully. No
#      trial spans both age groups, and each arm is defined by its trial. Two
#      arm labels do appear in both groups -- Peanut protein and Control, from
#      the two PPS strata -- and those genuinely are the same intervention in
#      both. They carry seven participants between them, one of them a single
#      older and a single younger control.
#   2. Does the design matrix say so? Only if the arm is labelled with its
#      trial. Pooled, the two PPS labels bridge the age groups and the model
#      fits at full rank while the age effect loads onto those coefficients.
#      This is the trap; the pooled fit is not safer for being estimable.
#   3. What does forcing it in cost? Aging stays estimable, so the cost is
#      precision rather than identifiability: the added columns inflate the
#      variance of the Aging contrast and spend residual degrees of freedom,
#      which moves a p-value distribution dense near the BH boundary. The
#      script measures the inflation rather than asserting a mechanism.
#   4. Where is supplement actually estimable? Inside one trial. NORE
#      randomised beetroot against placebo across 11 older adults, and PPS ran
#      peanut protein against control in each age stratum at n = 4 and n = 2.
#      NORE is tested here because it is the one two-arm randomisation with a
#      pre/post design in a single age group.
#
# Both codings are built from `supplement_cohort`, which collapses the three
# EAA arms into one cohort, "Ruple et al." That trial dosed whey and essential
# amino acids at roughly a tenth of what elicits a physiological response,
# restricted no one's dairy intake, and its food logs showed protein and EAA
# intake equal across arms. The randomisation exists on paper; it never
# delivered three distinct exposures, so splitting it three ways would spend
# two design columns on a contrast that was never administered.
#
# Collapsing it also removes the younger Placebo, which was never the same
# preparation as NORE's nitrate-depleted beetroot juice and shared no
# participant with it. What bridges the age groups afterwards is only the two
# PPS labels, which is the honest version of the picture.
#
# limma has one random-effect slot, spent on subject via duplicateCorrelation.
# There is no second slot for supplement; `variancePartition::dream` would be
# the tool if one were needed, but it cannot manufacture the missing overlap.

withr::local_dir(here::here())

pacman::p_load(
  readxl, readr, dplyr, tibble, stringr, tidyr, limma, proteoDA, openxlsx
)
source("04_Figures/shared/supplement_overview.R")
source("04_Figures/shared/devices.R")

set.seed(42)

# The node text panel A prints. The parent_study codes are what the metadata
# and every CSV record, and what a reader cross-references against S10 Table;
# they read as jargon on a figure, so the panel prints these and the workbook
# carries both. Each label names the parent trial because the claim the panel
# has to support is a trial-level one -- Beetroot's placebo and its juice read
# as unrelated exposures otherwise.
NODE_LABELS <- c(
  "Ruple et al." = "Ruple : pooled",
  "Beetroot juice" = "Beetroot : juice",
  "Placebo" = "Beetroot : placebo",
  "Peanut protein" = "Peanut : protein",
  "Control" = "Peanut : control",
  "Control/Placebo" = "Ctrl/Placebo"
)

OUT <- "03_DEP/c_data"
XLSX <- file.path(OUT, "03_DEP_results.xlsx")
RPT <- "03_DEP/b_reports/02_supplement"
FIG_PDF <- file.path(RPT, "pdf")
FIG_PNG <- file.path(RPT, "png")
for (d in c(RPT, FIG_PDF, FIG_PNG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

dal_norm <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
keep <- colnames(dal_norm$data)

meta <- read_excel("00_input/YvO_meta.xlsx") |>
  filter(Col_ID %in% keep) |>
  mutate(
    subject = str_remove(Col_ID, "_(Pre|Post)$"),
    group = factor(Group_Time,
      levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
    ),
    arm_full = factor(supplement),
    arm_pooled = factor(supplement_cohort),
    arm_ctrl = factor(if_else(
      supplement_cohort %in% c("Control", "Placebo"),
      "Control/Placebo", supplement_cohort
    )),
    arm_trial = factor(paste(parent_study, supplement_cohort, sep = ":"))
  )
meta <- meta[match(keep, meta$Col_ID), ]

message("Cohort by age group:")
print(table(meta$supplement_cohort, meta$Group))

# One participant per row, labelled with its trial, because the pooled label is
# what hides the nesting from the design matrix and the two are worth seeing
# side by side. This object is both the workbook cross-tab and panel A, so the
# sheet and the figure cannot drift apart.
by_age <- meta |>
  distinct(subject, parent_study, supplement, supplement_cohort, Group) |>
  count(parent_study, supplement, supplement_cohort, Group) |>
  complete(nesting(parent_study, supplement, supplement_cohort), Group,
    fill = list(n = 0L)
  ) |>
  pivot_wider(names_from = Group, values_from = n, values_fill = 0L) |>
  rename(recorded_arm = supplement, pooled_cohort = supplement_cohort) |>
  mutate(control_collapsed = if_else(
    pooled_cohort %in% c("Control", "Placebo"),
    "Control/Placebo", pooled_cohort
  )) |>
  mutate(total = Young + Old) |>
  arrange(desc(Young > 0), parent_study, recorded_arm) |>
  mutate(panel_label = unname(NODE_LABELS[pooled_cohort]))

message("\nCohort against age group, one row per recorded arm:")
print(as.data.frame(by_age))

# Question 1 and 2: what each coding tells the model.
rank_check <- function(term, label) {
  group <- meta$group
  mm <- model.matrix(~ 0 + group + term)
  dropped <- nonEstimable(mm)
  tibble(
    coding = label,
    columns = ncol(mm),
    rank = qr(mm)$rank,
    non_estimable = if (is.null(dropped)) "none" else paste(dropped, collapse = "; "),
    verdict = if (is.null(dropped)) {
      "fits silently despite the confound"
    } else {
      "limma drops a column"
    }
  )
}

design_diag <- bind_rows(
  rank_check(meta$arm_pooled, "supplement pooled across trials"),
  rank_check(meta$arm_full, "supplement as recorded across all arms"),
  rank_check(meta$arm_ctrl, "supplement with control and placebo collapsed"),
  rank_check(meta$arm_trial, "supplement labelled with its trial")
)
print(as.data.frame(design_diag))

# Question 3: fit the published model, then force supplement in under both
# codings, and compare. All recorded arms are kept alongside the pooled
# cohort because the submitted manuscript quotes them, and a reviewer should be
# able to see that the answer does not depend on which coding we chose.
fit_design <- function(design) {
  colnames(design) <- make.names(colnames(design))
  block_cor <- duplicateCorrelation(
    dal_norm$data, design,
    block = meta$subject
  )$consensus
  fit <- lmFit(dal_norm$data, design,
    block = meta$subject, correlation = block_cor
  )
  cm <- makeContrasts(
    Training_Young = groupYoung_Post - groupYoung_Pre,
    Training_Old = groupOld_Post - groupOld_Pre,
    Aging = groupOld_Pre - groupYoung_Pre,
    Interaction = (groupOld_Post - groupOld_Pre) - (groupYoung_Post - groupYoung_Pre),
    levels = design
  )
  eb <- eBayes(contrasts.fit(fit, cm), robust = TRUE)
  # A protein with too few observations leaves an NA coefficient, and
  # p.adjust propagates it.
  counts <- vapply(
    colnames(cm),
    \(cn) sum(p.adjust(eb$p.value[, cn], "BH") < 0.05, na.rm = TRUE),
    integer(1)
  )
  list(
    counts = counts,
    coef = eb$coefficients,
    p = eb$p.value,
    se = sqrt(eb$s2.post) * eb$stdev.unscaled,
    columns = ncol(design),
    rank = qr(design)$rank,
    df_residual = eb$df.residual[1]
  )
}

FITS <- c(
  published = "published (age collapsed)",
  arms = "supplement (all)",
  cohort = "supplement (Ruple pooled)",
  ctrl = "supplement (control collapsed)"
)

fits <- list(
  published = fit_design(model.matrix(~ 0 + group, meta)),
  arms = fit_design(model.matrix(~ 0 + group + arm_full, meta)),
  cohort = fit_design(model.matrix(~ 0 + group + arm_pooled, meta)),
  ctrl = fit_design(model.matrix(~ 0 + group + arm_ctrl, meta))
)

contrasts_in_order <- names(fits$published$counts)

forced <- tibble(
  contrast = contrasts_in_order,
  published = fits$published$counts,
  supplement_forced = fits$cohort$counts,
  supplement_forced_six_arms = fits$arms$counts,
  supplement_forced_control_collapsed = fits$ctrl$counts
) |>
  mutate(change = supplement_forced - published)

message("\nForcing supplement into the published model, both codings:")
print(as.data.frame(forced))

# The response letter quotes how far the coefficients move and how much the
# uncertainty grows, not just how many survive. Emit both here so the letter
# cites something the pipeline reproduces.
compare_to_published <- function(key) {
  lapply(contrasts_in_order, function(cn) {
    a <- fits$published$coef[, cn]
    b <- fits[[key]]$coef[, cn]
    sa <- fits$published$se[, cn]
    sb <- fits[[key]]$se[, cn]
    ok <- !is.na(a) & !is.na(b)
    tibble(
      contrast = cn,
      coding = unname(FITS[key]),
      n_compared = sum(ok),
      spearman_logfc = round(cor(a[ok], b[ok], method = "spearman"), 3),
      max_abs_shift = round(max(abs(a[ok] - b[ok])), 3),
      median_se_published = round(median(sa[ok]), 4),
      median_se_adjusted = round(median(sb[ok]), 4),
      se_ratio = round(median(sb[ok] / sa[ok]), 3)
    )
  }) |> bind_rows()
}

shift <- bind_rows(lapply(c("cohort", "arms", "ctrl"), compare_to_published))

message("\nCoefficient movement and uncertainty when supplement is forced in:")
print(as.data.frame(shift))

# Question 4: the NORE randomisation, one of three within-trial comparisons.
nore <- meta |> filter(parent_study == "NORE")
nore_mat <- dal_norm$data[, nore$Col_ID, drop = FALSE]

# Protein filtering was done across all 62 samples, so require presence in
# most of both arms before testing anything in this subset.
brj <- nore$supplement == "Beetroot juice"
stopifnot("both NORE arms present" = sum(brj) > 0 && sum(!brj) > 0)
# Two thirds of each arm, rather than counts tied to this cohort's sizes.
keep_p <- rowSums(!is.na(nore_mat[, brj])) >= ceiling(sum(brj) * 2 / 3) &
  rowSums(!is.na(nore_mat[, !brj])) >= ceiling(sum(!brj) * 2 / 3)
nore_mat <- nore_mat[keep_p, ]

nore_meta <- nore |>
  mutate(arm_time = factor(
    paste(if_else(brj, "BRJ", "PLA"), Timepoint, sep = "_"),
    levels = c("PLA_Pre", "PLA_Post", "BRJ_Pre", "BRJ_Post")
  )) |>
  as.data.frame()
rownames(nore_meta) <- nore_meta$Col_ID

dal <- DAList(
  data = nore_mat,
  annotation = dal_norm$annotation[keep_p, , drop = FALSE],
  metadata = nore_meta,
  tags = list(normalized = TRUE, norm_method = "cycloess")
) |>
  add_design("~ 0 + arm_time + (1 | subject)")

colnames(dal$design$design_matrix) <- gsub(
  "^arm_time", "", colnames(dal$design$design_matrix)
)

dal <- dal |>
  add_contrasts(contrasts_vector = c(
    "Baseline_arm   = BRJ_Pre - PLA_Pre",
    "Training_BRJ   = BRJ_Post - BRJ_Pre",
    "Training_PLA   = PLA_Post - PLA_Pre",
    "Arm_x_Training = (BRJ_Post - BRJ_Pre) - (PLA_Post - PLA_Pre)"
  )) |>
  fit_limma_model()

res <- extract_DA_results(dal, pval_thresh = 0.05, lfc_thresh = 0, adj_method = "BH")

nore_summary <- bind_rows(lapply(names(res$results), function(nm) {
  x <- res$results[[nm]]
  tibble(
    contrast = nm,
    n_tested = nrow(x),
    raw_p_lt_05 = sum(x$P.Value < 0.05, na.rm = TRUE),
    expected_by_chance = round(0.05 * nrow(x)),
    fdr_lt_05 = sum(x$adj.P.Val < 0.05, na.rm = TRUE),
    min_adj_p = signif(min(x$adj.P.Val, na.rm = TRUE), 3)
  )
}))

message("\nNORE beetroot against placebo, randomised within the older cohort:")
print(as.data.frame(nore_summary))

write_csv(design_diag, file.path(OUT, "02_supplement_design_diagnostics.csv"))
write_csv(forced, file.path(OUT, "02_supplement_forced_fit.csv"))
write_csv(shift, file.path(OUT, "02_supplement_forced_shift.csv"))
write_csv(nore_summary, file.path(OUT, "02_supplement_nore_randomised.csv"))

# The workbook is a table artifact, so it lands before the reports for the same
# reason the CSVs do. This is the last script to write 03_DEP_results.xlsx, so
# it owns the Overview index.


wb <- loadWorkbook(XLSX)
write_sheet(wb, "supplement_by_age", by_age)
write_sheet(wb, "supplement_design", design_diag)
write_sheet(wb, "supplement_forced_fit", forced)
write_sheet(wb, "supplement_forced_shift", shift)
write_sheet(wb, "supplement_nore", nore_summary)

add_overview(
  wb,
  title = "S10 Table — differential abundance",
  description = paste(
    "Every protein in every contrast, the robustness checks behind the counts",
    "the paper reports, and the supplementation analysis answering R1.2."
  ),
  entries = c(
    combined_results        = "One row per protein with the fold change, p-value, adjusted p-value and Π-score for all four contrasts.",
    DA_summary              = "How many proteins reached each significance threshold, per contrast and direction.",
    Aging                   = "Older against younger at baseline, ranked by Π-score.",
    Training_Young          = "Post against pre in younger adults, ranked by Π-score.",
    Training_Old            = "Post against pre in older adults, ranked by Π-score.",
    Interaction             = "Whether the training response differs by age group, ranked by Π-score.",
    outlier_sensitivity     = "Whether the counts hold when the two outlier samples are kept in.",
    blunting                = "How far the older-adult effect sizes are compressed against the younger.",
    bootstrap_ci            = "Confidence intervals on the headline counts, by resampling participants.",
    power_analysis          = "Post-hoc power for each contrast at the observed effect sizes.",
    imputation_sensitivity  = "Whether the counts hold when the imputation method is changed.",
    supplement_by_age       = "Participants per recorded arm per age group, with the cohort label each arm carries once the three EAA arms are pooled, and the node label panel A of S9 Figure prints for that cohort. No trial appears in both age groups; only peanut protein and control bridge them, and each control cell holds one participant.",
    supplement_design       = "Design matrix columns, rank and non-estimable columns for three codings of the supplement term: pooled to parent-trial cohort, as recorded across all arms, with control and placebo collapsed, and labelled with its trial.",
    supplement_forced_fit   = "Proteins at FDR < 0.05 per contrast for the published model and for the same model with supplement forced in, under all three adjusted codings: all recorded arms, Ruple pooled, and control and placebo collapsed.",
    supplement_forced_shift = "Spearman correlation, largest absolute shift and median standard error between the published and forced fits, per contrast and per coding.",
    supplement_nore         = "Beetroot juice against placebo within NORE, the one two-arm randomisation with a pre/post design in a single age group."
  )
)
saveWorkbook(wb, XLSX, overwrite = TRUE)

# Reports last: rendering is the fragile half and must not cost the tables.
write_limma_plots(res,
  grouping_column = "arm_time", output_dir = RPT, overwrite = TRUE
)
write_limma_tables(res, output_dir = RPT, overwrite = TRUE)

message("\nDone -> ", OUT, " and ", RPT)

# Four panels: who was in what, what adjusting for it costs, what it does to
# the uncertainty, and whether the effect sizes themselves move. Both forced
# codings are drawn throughout, because the conclusion should not rest on which
# one we picked. Every panel reads an object built above, so the figure cannot
# drift from the tables it ships with.
pacman::p_load(ggplot2, patchwork)

BLUE <- "#4393C3"
SAND <- "#F4A582"
RED <- "#D6604D"
GREY <- "grey35"


# Short forms for panel A's column heads and panel D's row strips. The full
# names are in panel B's legend; at their length they clip inside a strip and
# collide across the column heads.
SHORT_NAMES <- c(
  published = "age collapsed",
  arms = "all arms",
  cohort = "Ruple pooled",
  ctrl = "ctrl/placebo"
)

# The all-arms coding is fitted and tabled but never drawn. Pooling by
# supplement label across trials fuses NORE's nitrate-depleted placebo with
# Ruple's into one ten-person level, and setting that partition beside the
# others invited the reader to weigh it as a candidate design. Its counts stay
# in supplement_forced_fit.
SHOWN <- c("published", "cohort", "ctrl")

STRIP_LABELS <- setNames(
  unname(SHORT_NAMES[c("cohort", "ctrl")]),
  FITS[c("cohort", "ctrl")]
)

# Blue is the published reference; the two adjusted fits run light to dark in
# the order they coarsen, so the ramp itself reads as "more pooling".
FIT_COLOURS <- c(
  `published (age collapsed)` = BLUE,
  `supplement (Ruple pooled)` = RED,
  `supplement (control collapsed)` = "#7F2A17"
)

# Drawn at the width it is read at, so declared point sizes are printed point
# sizes and tools/check_figures.py measures them at a scale of 1.
FIG_W <- 178
FIG_H <- 275

base_theme <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    legend.position = "none"
  )

contrast_label <- function(x) str_replace(x, "_", " in ")

# Panel A reads the four fitted models as successive coarsenings of one
# partition. The leftmost column is the data as recorded, trial by arm, and is
# not a model: "all arms" pools by supplement label across trials, which is
# where Placebo from two different trials becomes one ten-person level that

# Panel A reads the three fitted models as successive coarsenings of one
# partition. The leftmost column is the coarsening the robustness check starts
# from: Ruple's three arms pooled, every other arm as recorded.
stopifnot(
  "every pooled cohort has a panel A label" =
    all(by_age$pooled_cohort %in% names(NODE_LABELS)),
  "every collapsed level has a panel A label" =
    all(by_age$control_collapsed %in% names(NODE_LABELS))
)
short_node <- function(x) if_else(x %in% names(NODE_LABELS), NODE_LABELS[x], x)

# One row per cohort per age group. A cohort spanning both groups needs the
# split to draw its edges, while the node itself is drawn once.
base <- by_age |>
  select(pooled_cohort, control_collapsed, Young, Old) |>
  pivot_longer(c(Young, Old), names_to = "Group", values_to = "n") |>
  filter(n > 0) |>
  count(pooled_cohort, control_collapsed, Group, wt = n, name = "n")

SPAN <- c(1, 5)

# Spaced evenly rather than at the mean of their children: a node drawn from
# both age groups lands on the same mid-point as one drawn from a single
# group, and the labels collide. Ordered young-only, then the ones spanning
# both ages, then old-only, so the streams that cross sit in the middle, and
# within a band by where the children sit, which is what keeps the streams
# from crossing at all.
spread <- function(d) {
  d |>
    mutate(age_rank = match(span, c("Young", "both", "Old"))) |>
    arrange(age_rank, desc(ord)) |>
    mutate(y = seq(max(SPAN), min(SPAN), length.out = n()))
}

lvl1 <- base |>
  group_by(node = pooled_cohort) |>
  summarise(
    ord = sum(n),
    n = sum(n),
    both_ages = n_distinct(Group) > 1,
    span = if (n_distinct(Group) > 1) "both" else Group[1],
    .groups = "drop"
  ) |>
  spread()

collapse_level <- function(key) {
  base |>
    left_join(select(lvl1, pooled_cohort = node, child_y = y),
      by = "pooled_cohort"
    ) |>
    group_by(node = .data[[key]]) |>
    summarise(
      ord = mean(child_y),
      n = sum(n),
      both_ages = n_distinct(Group) > 1,
      span = if (n_distinct(Group) > 1) "both" else Group[1],
      .groups = "drop"
    ) |>
    spread()
}

lvl2 <- collapse_level("control_collapsed")
lvl3 <- collapse_level("Group")

step_edges <- function(from_key, to_key, from_tbl, to_tbl) {
  base |>
    count(
      .from = .data[[from_key]], .to = .data[[to_key]], wt = n, name = "w"
    ) |>
    left_join(select(from_tbl, node, y_from = y), by = c(.from = "node")) |>
    left_join(select(to_tbl, node, y_to = y, both_ages), by = c(.to = "node"))
}

edge12 <- step_edges("pooled_cohort", "control_collapsed", lvl1, lvl2)
edge23 <- step_edges("control_collapsed", "Group", lvl2, lvl3)

node_label <- function(label, n) sprintf("%s  (%d)", short_node(label), n)

FORMULAE <- c(
  cohort = "~ 0 + group + arm(5)",
  ctrl = "~ 0 + group + arm(4)",
  published = "~ 0 + group"
)

design_header <- function(key) {
  f <- fits[[key]]
  sprintf(
    "%s\n%s\n%d cols, rank %d\n%d resid df",
    SHORT_NAMES[[key]], FORMULAE[[key]], f$columns, f$rank, f$df_residual
  )
}

step_layer <- function(d, x) {
  geom_segment(
    data = d,
    aes(
      x = x, xend = x + 1, y = y_from, yend = y_to,
      colour = both_ages, linewidth = w
    )
  )
}

node_layer <- function(d, x, hjust = 0.5, bold = TRUE) {
  geom_label(
    data = d,
    aes(x = x, y = y, label = node_label(node, n), colour = both_ages),
    hjust = hjust, size = 2.2, label.size = 0.2,
    label.padding = unit(0.7, "mm"), fill = "white",
    fontface = if (bold) "bold" else "plain"
  )
}

p_design <- ggplot() +
  step_layer(edge12, 1) +
  step_layer(edge23, 2) +
  node_layer(lvl1, 1, hjust = 1, bold = FALSE) +
  node_layer(lvl2, 2) +
  node_layer(lvl3, 3, hjust = 0) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = RED)) +
  scale_linewidth_continuous(range = c(0.25, 1.5)) +
  scale_x_continuous(
    # The left gutter has to hold the longest cohort label, which is wider
    # than the gap between node columns; too narrow and ggplot clips the first
    # character inside the panel, where a span check cannot see it.
    limits = c(0.45, 3.5),
    breaks = 1:3,
    labels = vapply(
      c("cohort", "ctrl", "published"), design_header, character(1)
    ),
    position = "top"
  ) +
  labs(
    title = "Cohort structure and design under each model",
    x = NULL, y = NULL,
    caption = paste(
      "cols, parameters the model estimates; rank, how many carry independent",
      "information. Equal means every term is estimable.\nRed marks a level",
      "holding participants from both age groups."
    )
  ) +
  base_theme +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      size = 6.6, face = "bold", colour = "grey25",
      lineheight = 1.2
    ),
    plot.caption = element_text(
      size = 6.2, colour = "grey40", hjust = 0,
      lineheight = 1.3
    )
  )


counts_long <- forced |>
  select(contrast,
    `published (age collapsed)` = published,
    `supplement (Ruple pooled)` = supplement_forced,
    `supplement (control collapsed)` = supplement_forced_control_collapsed
  ) |>
  pivot_longer(-contrast, names_to = "fit", values_to = "n") |>
  mutate(
    contrast = factor(contrast_label(contrast),
      levels = contrast_label(contrasts_in_order)
    ),
    fit = factor(fit, levels = names(FIT_COLOURS))
  )

p_counts <- ggplot(counts_long, aes(contrast, n, fill = fit)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.78) +
  geom_text(aes(label = n),
    position = position_dodge(width = 0.85), vjust = -0.35, size = 2.05
  ) +
  scale_fill_manual(values = FIT_COLOURS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(
    title = "Proteins at FDR < 0.05",
    x = NULL, y = "proteins", fill = NULL
  ) +
  base_theme +
  theme(
    legend.position = "top",
    legend.key.size = unit(3, "mm"),
    legend.text = element_text(size = 7),
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

se_long <- bind_rows(lapply(SHOWN, function(k) {
  se <- fits[[k]]$se
  tibble(
    fit = unname(FITS[k]),
    contrast = rep(colnames(se), each = nrow(se)),
    se = as.vector(se)
  )
})) |>
  filter(!is.na(se)) |>
  mutate(
    contrast = factor(contrast_label(contrast),
      levels = contrast_label(contrasts_in_order)
    ),
    fit = factor(fit, levels = names(FIT_COLOURS))
  )

# Outliers are suppressed and the axis clipped to the boxes: the quartiles are
# computed on every protein, and a long tail of high-variance proteins common
# to all three fits would otherwise flatten the difference that matters.
p_se <- ggplot(se_long, aes(contrast, se, fill = fit)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
  scale_fill_manual(values = FIT_COLOURS) +
  coord_cartesian(ylim = c(0, 0.45)) +
  labs(
    title = "Standard error by contrast",
    x = NULL, y = "standard error"
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

# One point per protein, published against adjusted. Points on the identity
# line did not move. The three within-age contrasts are the control: they sit
# on the line under both codings, so the aging scatter is not a general effect
# of adding columns.
scatter_long <- bind_rows(lapply(c("cohort", "ctrl"), function(k) {
  tibble(
    fit = unname(FITS[k]),
    contrast = rep(contrasts_in_order, each = nrow(fits$published$coef)),
    published = as.vector(fits$published$coef[, contrasts_in_order]),
    adjusted = as.vector(fits[[k]]$coef[, contrasts_in_order])
  )
})) |>
  filter(!is.na(published), !is.na(adjusted)) |>
  mutate(
    contrast = factor(contrast_label(contrast),
      levels = contrast_label(contrasts_in_order)
    ),
    fit = factor(fit, levels = names(FIT_COLOURS))
  )

rho_labels <- shift |>
  filter(coding %in% FITS[c("cohort", "ctrl")]) |>
  transmute(
    contrast = factor(contrast_label(contrast),
      levels = contrast_label(contrasts_in_order)
    ),
    fit = factor(coding, levels = names(FIT_COLOURS)),
    label = sprintf("rho == %.3f", spearman_logfc)
  )

p_scatter <- ggplot(scatter_long, aes(published, adjusted, colour = fit)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.25, alpha = 0.3) +
  geom_text(
    data = rho_labels, parse = TRUE,
    aes(x = -Inf, y = Inf, label = label),
    hjust = -0.12, vjust = 1.4, size = 2.3, show.legend = FALSE
  ) +
  # Short strip labels: the full model names are already in panel B's legend,
  # and at 25 characters the row strip clips them inside the panel, where the
  # span check cannot see it.
  facet_grid(fit ~ contrast, labeller = labeller(fit = STRIP_LABELS)) +
  scale_x_continuous(breaks = c(-2, 0, 2)) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  scale_colour_manual(values = FIT_COLOURS) +
  labs(
    title = "Adjusted against published fold change",
    x = "published log fold change", y = "adjusted log fold change"
  ) +
  base_theme +
  theme(
    strip.text = element_text(size = 7),
    panel.spacing = unit(2, "mm")
  )

# Panel E: where the counts come from. A real effect piles p-values against
# zero and a null spreads them flat, so the aging signal draining away reads as
# a shape rather than as four numbers that could be dismissed as an artefact of
# where the threshold sits.
pval_long <- bind_rows(lapply(SHOWN, function(k) {
  tibble(fit = unname(FITS[k]), p = fits[[k]]$p[, "Aging"])
})) |>
  filter(!is.na(p)) |>
  mutate(fit = factor(fit, levels = names(FIT_COLOURS)))

p_pvals <- ggplot(pval_long, aes(p, fill = fit)) +
  geom_histogram(
    binwidth = 0.025, boundary = 0, colour = "white",
    linewidth = 0.1
  ) +
  geom_vline(
    xintercept = 0.05, colour = GREY, linetype = "22",
    linewidth = 0.3
  ) +
  facet_wrap(~fit, nrow = 1) +
  scale_fill_manual(values = FIT_COLOURS) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  labs(
    title = "Aging p-values: the signal flattens as the design coarsens",
    x = "raw p-value", y = "proteins"
  ) +
  base_theme +
  theme(strip.text = element_text(size = 6.8), panel.spacing = unit(2, "mm"))

fig <- (p_design / (p_counts | p_se) / p_scatter / p_pvals) +
  plot_layout(heights = c(0.95, 0.95, 0.95, 0.7)) +
  plot_annotation(
    tag_levels = "A",
    title = paste(
      "No parent trial spans both age groups, so adjusting for cohort",
      "removes the aging signal"
    ),
    theme = theme(plot.title = element_text(face = "bold", size = 10))
  )

# open_pdf resolves the device rather than letting ggsave pick base pdf from
# the extension, which drops characters the other staged figures keep.
open_pdf(
  file.path(FIG_PDF, "SUPP_DEP_covariate.pdf"), FIG_W / 25.4, FIG_H / 25.4
)
print(fig)
dev.off()

ggsave(file.path(FIG_PNG, "SUPP_DEP_covariate.png"), fig,
  width = FIG_W, height = FIG_H, units = "mm", dpi = 300, bg = "white"
)

message("\nFigure -> ", FIG_PDF, " and ", FIG_PNG)
