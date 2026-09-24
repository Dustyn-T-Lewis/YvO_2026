#!/usr/bin/env Rscript
# Figure 6C: module x phenotype coupling is age-dependent (3x2 grid).
#
# Six hero scatters, taken as the six best source x module x outcome cells in
# the live screen; each cell plots both age strata. Writes the 180-test screen
# to c_data/panel_B_full_screen_bh.csv for S7_B_full_sweep.R and S7 Table.

setwd(here::here())

pacman::p_load(tidyverse, patchwork, ggtext)

source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")

BASE <- "04_Figures/F06"
RPT <- file.path(BASE, "b_reports", "panels")
DAT_OUT <- file.path(BASE, "c_data")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_OUT, recursive = TRUE, showWarnings = FALSE)

F05_SUPP <- "04_Figures/F05/c_data/F05_data.xlsx"
stopifnot(
  "run 04_Figures/F05/a_script/F05_data.R first: missing F05_data.xlsx" =
    file.exists(F05_SUPP)
)

me_pre <- read_matrix_sheet(F05_SUPP, "me_pre", "subject_key")
delta_me <- read_matrix_sheet(F05_SUPP, "delta_me", "subject_key")
pheno_wide <- read_sheet_df(F05_SUPP, "metadata_pheno_wide")
subj_age <- read_sheet_df(F05_SUPP, "metadata_subj_age")
common_subj <- read_vector_sheet(F05_SUPP, "common_subj")

pdf_device <- get_pdf_device()

message("Panel B: 3x2 grid, top 6 cells of the module x phenotype screen...")

# Outcome + predictor labels
# Spelled out because T1/T2/DL are storage names, not things a reader knows:
# T1/T2 are Type I / Type II fibre cross-sectional area, DL is deadlift 1RM.
outcome_nice <- c(
  delta_VL = "\u0394VL thickness (cm)",
  delta_LBM = "\u0394Lean mass (kg)",
  delta_DL = "\u0394Deadlift 1RM (kg)",
  delta_T1 = "\u0394Type I fCSA (\u00b5m\u00b2)",
  delta_T2 = "\u0394Type II fCSA (\u00b5m\u00b2)"
)

outcome_short <- c(
  delta_VL = "\u0394VL",
  delta_LBM = "\u0394LBM",
  delta_DL = "\u0394Deadlift",
  delta_T1 = "\u0394Type I fCSA",
  delta_T2 = "\u0394Type II fCSA"
)

# The colour name, not the pathway. A module is named by its colour, the cell
# is framed in that colour, and the figure key carries colour -> pathway once
# so six cells do not each spell it out.
pretty_mod <- function(m) str_to_title(gsub("^ME", "", m))

get_x <- function(x_source, mod) {
  src <- if (x_source == "baseline") me_pre else delta_me
  as.numeric(src[common_subj, mod])
}

# Helpers
fmt_p <- function(p) {
  if (is.na(p)) {
    return("n/a")
  }
  if (p < 0.001) "p<0.001" else sprintf("p=%.3f", p)
}
sig_mark <- function(p) {
  if (is.na(p)) {
    ""
  } else if (p < 0.001) {
    "***"
  } else if (p < 0.01) {
    "**"
  } else if (p < 0.05) {
    "*"
  } else if (p < 0.10) {
    "\u2020"
  } else {
    ""
  }
}

# Full screening + BH correction audit
MODULES_SCREEN <- grep("^ME", colnames(delta_me), value = TRUE)
MODULES_SCREEN <- setdiff(MODULES_SCREEN, "MEgrey")
SOURCES <- c("delta_ME", "baseline")
OUTCOMES <- names(outcome_nice)

screen_rows <- list()
for (src in SOURCES) {
  src_mat <- if (src == "baseline") me_pre else delta_me
  for (mod in MODULES_SCREEN) {
    x <- as.numeric(src_mat[common_subj, mod])
    for (out in OUTCOMES) {
      y <- pheno_wide[[out]][match(common_subj, pheno_wide$subject_key)]
      for (grp in c("Young", "Old")) {
        age <- subj_age$age[match(common_subj, subj_age$subject_key)]
        keep <- !is.na(x) & !is.na(y) & age == grp
        if (sum(keep) < 4) next
        ct <- suppressWarnings(cor.test(x[keep], y[keep]))
        screen_rows[[length(screen_rows) + 1]] <- tibble(
          source = src, module = mod, outcome = out, stratum = grp,
          n = sum(keep), r = unname(ct$estimate), p_raw = ct$p.value
        )
      }
    }
  }
}

screen_df <- bind_rows(screen_rows) |>
  mutate(p_bh = p.adjust(p_raw, method = "BH")) |>
  arrange(p_raw)

n_screen <- nrow(screen_df)
n_raw_sig <- sum(screen_df$p_raw < 0.05)
n_bh_sig <- sum(screen_df$p_bh < 0.05)
message(sprintf(
  "  Full screen: %d tests, %d raw p<0.05, %d BH p<0.05",
  n_screen, n_raw_sig, n_bh_sig
))

write_csv(screen_df, file.path(DAT_OUT, "panel_B_full_screen_bh.csv"))

N_HERO <- 6L

# The panel used to carry a hand-curated pick-list, which drifted: two of its
# four entries were not even nominally significant by the time this ran. Take
# the six strongest cells from the screen itself instead. A cell is one
# source x module x outcome and is drawn with both age strata, so it is ranked
# on whichever stratum does better.
HERO_PICKS <- screen_df |>
  group_by(source, module, outcome) |>
  slice_min(p_raw, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(p_raw) |>
  slice_head(n = N_HERO) |>
  transmute(x_source = source, module, outcome)

# The old version of this check compared nrow() of a stratum-expanded join
# against the number of picks. Each surviving pick contributed a Young row and
# an Old row, so two passing picks masked two failing ones and the warning
# never fired. Count distinct cells, not joined rows.
hero_keys <- with(HERO_PICKS, paste(x_source, module, outcome, sep = "|"))
passing_keys <- screen_df |>
  filter(p_raw < 0.10) |>
  transmute(key = paste(source, module, outcome, sep = "|")) |>
  pull(key) |>
  unique()
failed_keys <- setdiff(hero_keys, passing_keys)
if (length(failed_keys)) {
  warning(sprintf(
    "%d of %d hero picks have no stratum with raw p<0.10: %s",
    length(failed_keys), length(hero_keys),
    paste(failed_keys, collapse = ", ")
  ))
}

# Stats come from screen_df rather than a second cor.test so the cell text and
# the exported screen can never disagree, and so q is available at all.
screen_stat <- function(x_source, mod, out, grp) {
  hit <- screen_df[screen_df$source == x_source & screen_df$module == mod &
    screen_df$outcome == out & screen_df$stratum == grp, ]
  if (nrow(hit) != 1L) {
    return(c(r = NA_real_, p = NA_real_, q = NA_real_, n = NA_real_))
  }
  c(r = hit$r[1], p = hit$p_raw[1], q = hit$p_bh[1], n = hit$n[1])
}

build_hero_mini <- function(x_source, module, outcome) {
  mod <- module
  x <- get_x(x_source, mod)
  y <- pheno_wide[[outcome]][match(common_subj, pheno_wide$subject_key)]
  age <- subj_age$age[match(common_subj, subj_age$subject_key)]

  d <- tibble(x = x, y = y, age = factor(age, levels = c("Young", "Old"))) |>
    filter(!is.na(x), !is.na(y))

  cY <- screen_stat(x_source, mod, outcome, "Young")
  cO <- screen_stat(x_source, mod, outcome, "Old")

  # q is on every line because every cell in this screen has q = 0.94: the
  # panel has to say out loud that none of it survives BH correction.
  stat_line <- function(tag, st) {
    if (is.na(st[["r"]])) {
      return(sprintf("%s: too few subjects", tag))
    }
    # The dagger collides with the comma that follows it unless spaced off.
    mark <- sig_mark(st[["p"]])
    if (nzchar(mark)) mark <- paste0("\u2009", mark)
    sprintf(
      "%s (n=%d): r=%+.2f, %s%s, q=%.2f",
      tag, as.integer(st[["n"]]), st[["r"]],
      fmt_p(st[["p"]]), mark, st[["q"]]
    )
  }
  y_stat <- stat_line("Y", cY)
  o_stat <- stat_line("O", cO)

  mod_hex <- module_fill(mod)
  # Pre and Pre-Post, the vocabulary panels A and B already use, rather than
  # "baseline" and a bare delta. The y axis carries its own delta, so the delta
  # here names what it is a change in.
  #
  # A hyphen, not an arrow: U+2192 survives cowplot::draw_label but drops out of
  # annotate() inside the panel, which subsets a font without it, and the glyph
  # vanishes with no warning.
  x_prefix <- if (x_source == "baseline") "Pre" else "\u0394 Pre-Post"
  x_lab <- "Eigengene value (a.u.)"
  y_lab <- unname(outcome_nice[outcome])
  # The module bio-labels run to 30 characters, so the cell title takes the
  # short outcome token; the spelled-out version stays on the y axis.
  title <- sprintf(
    "%s %s vs %s",
    x_prefix, pretty_mod(mod), unname(outcome_short[outcome])
  )

  stat_label <- sprintf(
    "<span style='color:%s'><b>%s</b></span><br><span style='color:%s'><b>%s</b></span>",
    AGE_COLORS[["Young"]], y_stat, AGE_COLORS[["Old"]], o_stat
  )

  ggplot(d, aes(x = x, y = y, color = age, fill = age)) +
    geom_hline(
      yintercept = 0, linetype = "dotted",
      color = "grey70", linewidth = 0.2
    ) +
    geom_smooth(
      data = filter(d, age == "Young"),
      method = "lm", se = TRUE, linewidth = 0.6, alpha = 0.13
    ) +
    geom_smooth(
      data = filter(d, age == "Old"),
      method = "lm", se = TRUE, linewidth = 0.6,
      alpha = 0.11, linetype = "dashed"
    ) +
    geom_point(
      shape = 21, size = 2.4, stroke = 0.4,
      color = "white", alpha = 0.95
    ) +
    scale_color_manual(values = AGE_COLORS) +
    scale_fill_manual(values = AGE_COLORS) +
    scale_y_continuous(expand = expansion(mult = c(0.25, 0.12))) +
    labs(x = x_lab, y = y_lab) +
    # Title inside plot area at top centre, in the module's own colour
    # darkened to stay readable as type. The name carries the module in text as
    # well as in the frame colour, which is the cue a colourblind reader loses.
    annotate("text",
      x = mean(range(d$x, na.rm = TRUE)), y = Inf,
      label = title, hjust = 0.5, vjust = 1.15,
      size = 2.4, fontface = "bold", color = text_safe(mod_hex)
    ) +
    ggtext::geom_richtext(
      inherit.aes = FALSE,
      data = tibble(x = mean(range(d$x, na.rm = TRUE)), y = -Inf, lab = stat_label),
      aes(x = x, y = y, label = lab),
      hjust = 0.5, vjust = -0.1,
      size = 2.2, lineheight = 1.35,
      label.color = "grey70", label.r = unit(1.2, "pt"),
      label.padding = unit(c(1.5, 2.5, 1.5, 2.5), "pt"),
      label.margin = unit(c(2, 2, 2, 2), "pt"),
      fill = scales::alpha("white", 0.88)
    ) +
    FIG_THEME +
    coord_cartesian(clip = "off") +
    theme(
      axis.title.x = element_text(
        size = 5, face = "bold",
        margin = margin(t = 3)
      ),
      # grid clamps a negative right margin here, so -6 and -10 render
      # byte-identically; the value is inert. The fCSA titles still graze
      # the -2000 tick by 8.3 pt-squared. Fixing that needs the tick itself
      # shortened, e.g. scales::label_number(scale = 1e-3, suffix = "k").
      axis.title.y = element_text(
        size = 5, face = "bold",
        margin = margin(r = -6, unit = "pt")
      ),
      axis.text = element_text(size = FIG_AXIS_TEXT),
      legend.position = "none",
      # Thin: the frame only has to say which module the cell belongs to, and a
      # heavy one competes with the regression lines inside it.
      panel.border = element_rect(colour = mod_hex, linewidth = 0.45, fill = NA),
      plot.margin = margin(4, 2, 4, 2)
    )
}

hero_plots <- pmap(HERO_PICKS, build_hero_mini)

GRID_NCOL <- 3L # 6 picks across a 210 mm canvas: 3 wide keeps cells legible

# Remove x-axis title from all but the bottom row
for (i in seq_along(hero_plots)) {
  if (i <= length(hero_plots) - GRID_NCOL) {
    hero_plots[[i]] <- hero_plots[[i]] +
      theme(axis.title.x = element_blank())
  }
}

# Naming the BH family matters: Figure 5 corrects the same raw p within a
# trait column (10 modules) and reports q=0.058 for the green/dVL Young cell,
# where this screen's 180-test family gives q=0.94.
pC_subtitle <- sprintf(
  paste(
    "Pre or \u0394 Pre-Post ME vs \u0394 phenotype | r by stratum |",
    "BH over %d tests: %d raw, %d q | labels: enrichment then hub"
  ),
  n_screen, n_raw_sig, n_bh_sig
)

panel_B <- wrap_plots(hero_plots, ncol = GRID_NCOL) +
  plot_annotation(
    title = "Module\u2013Phenotype Coupling (Age-Dependent)",
    subtitle = pC_subtitle,
    theme = theme(
      plot.title = element_text(
        face = "bold", size = FIG_TITLE_SIZE,
        lineheight = 1.15,
        margin = margin(
          t = 6, l = 34, b = 3,
          unit = "pt"
        )
      ),
      plot.subtitle = element_text(
        size = FIG_SUBTITLE_SIZE,
        face = "bold.italic", color = "grey30",
        lineheight = 1.2,
        margin = margin(
          t = 2, l = 34, b = 6,
          unit = "pt"
        )
      )
    )
  )

PC_W <- 210 # matches the width panel C gets in the composite
PC_H <- 88

ggsave(file.path(RPT, "C_hero_grid.png"),
  panel_B,
  width = PC_W, height = PC_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT, "C_hero_grid.pdf"),
  panel_B,
  width = PC_W, height = PC_H, units = "mm",
  device = get_pdf_device()
)

message("  C_hero_grid saved (3x2 hero grid)")

invisible(panel_B)
