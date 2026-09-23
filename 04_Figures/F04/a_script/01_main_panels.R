#!/usr/bin/env Rscript
# F04 Main: training plasticity across age (6-panel composite)
#
# A  protein quadrant scatter + flanking quadrant ORA, Tr.(Y) x Tr.(O)
# B  RRHO2 threshold-free rank overlap, same pair
# C  response magnitude and coherence per subject
# D  pathway NES scatter, Tr.(Y) x Tr.(O)
# E  pathway NES scatter, Aging x Tr.(O)
# F  this study's own FDR and Pi sets against the Training(Old) ranking
#
# A, B, C and D use contrast pairs sharing no group mean. E pairs Aging with
# Training(Old), which share Old_Pre; pathway aggregation shrinks that coupling
# because the shared per-protein noise averages down while genuine
# between-pathway differences do not, so E is reported at pathway level only.
# F's Aging rows are coupled at protein level and take a Pre/Post swap null.
#
# F05 was folded into this figure on 2026-08-24; see docs/decisions.md.
# limma tree: 2106 proteins over 62 samples, and Training(Old) has no FDR hits,
# so panel F's Tr.(O) rows are empty by construction rather than by filtering.

setwd(here::here())

pacman::p_load(
  dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot, readxl
)

# ggrepel places labels by a stochastic search, so an unseeded render puts
# them somewhere new each time. run_all.R runs each script in its own
# Rscript child, which starts from a time-seeded RNG.
set.seed(42)

source("04_Figures/shared/style.R")

BASE <- "04_Figures/F04"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT <- file.path(BASE, "c_data")
# Shared engines reassign RPT_PNG/RPT_PDF to their panels subdir, so the
# composite writes through names they never touch.
OUT_PNG <- RPT_PNG
OUT_PDF <- RPT_PDF
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, SUPP_PNG, SUPP_PDF, DAT)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

pdf_device <- get_pdf_device()

# Shared engines take loaded objects, not paths, so the limma tree can supply
# the same structures from its own DAList without the engines knowing which
# pipeline they are running under.
# limma drops outliers upstream, so the metadata sheet carries samples the
# imputed matrix no longer has. Subset to the intersection before anything
# downstream builds a design from it.
abundance <- readRDS("02_imputation/c_data/01_DAList_imputed.rds")$data
meta <- as.data.frame(read_excel("00_input/YvO_meta.xlsx"))
meta <- meta[meta$Col_ID %in% colnames(abundance), ]

engine_cfg <- list(
  meta = meta,
  matrix = as.matrix(abundance[, meta$Col_ID]),
  dep_df = read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE),
  group_time = factor(meta$Group_Time,
    levels = c("Young_Pre", "Young_Post", "Old_Pre", "Old_Post")
  )
)

PATHWAY_LABELS <- c(
  "Unfolded Protein Response" = "UPR",
  "Ribosome Biogenesis" = "Ribo Bio",
  "Amino Acid Metabolism" = "AA Metabolism",
  "Amino Acid Metabolic Process" = "AA Metabolism",
  "Modified Amino Acid Metabolic Process" = "Modified AA Metabolism",
  "Sulfur Compound Metabolic Process" = "Sulfur Cmpd Metabolism",
  "Fatty Acid Metabolism" = "FA Metabolism",
  "Fatty Acid Beta Oxidation" = "FA Beta-Oxidation",
  "Oxidative Phosphorylation" = "OXPHOS",
  "Epithelial Mesenchymal Transition" = "EMT",
  "Plasma Membrane Protein Loc." = "PM Protein Loc.",
  "Cytoplasmic Translation" = "Cytoplasmic Transl.",
  "Mitochondrial Organization" = "Mito Org.",
  "Mitochondrion Organization" = "Mito Organization",
  "Mitochondrial Gene Expression" = "Mito Gene Expression",
  "Generation Of Precursor Metabolites And Energy" = "Precursor Metab. & Energy",
  "Precursor Metabolites & Energy" = "Precursor Metab. & Energy",
  "Extracellular Matrix Organization" = "ECM Organization",
  "Microtubule-Based Movement" = "MT-Based Movement",
  "Trna Metabolic Process" = "tRNA Metabolism",
  "Muscle System Process" = "Muscle System",
  "Cholesterol Homeostasis" = "Cholesterol Homeost.",
  "Heme Metabolism" = "Heme Metab.",
  "Ketone Metabolism" = "Ketone Metab."
)

message("F04 composite: sourcing panels")

source("04_Figures/F04/a_script/_panel_A_ORA.R")
panel_a <- composite
n_total_A <- nrow(scatter_df)
n_sig_A <- n_sig

cfg <- list(
  fig_id = "F04",
  t_col_1 = "t_Training_Young",
  t_col_2 = "t_Training_Old",
  rrho_labels = c("Training (Young)", "Training (Old)"),
  title = "Threshold-Free Concordance (RRHO2)",
  subtitle_fmt = "%s | %d shared proteins | warm corners = concordant regulation",
  axis_label_1 = expression("Training (Young) rank" ~ (Down %->% Up)),
  axis_label_2 = expression("Training (Old) rank" ~ (Down %->% Up)),
  # Position says which discordance is which, the way panels D and E already
  # label both their off-diagonal corners.
  quadrant_labels = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant", DU = "Discordant"
  ),
  hotspot_export_names = list(UU = "UU", DD = "DD", UD = "UD", DU = "DU"),
  ora_min_size = 15,
  ora_quadrant_names = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (Y Up / O Down)", DU = "Discordant (Y Down / O Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UU", "ora_DD"),
    file_2_quads = c("ora_UD", "ora_DU"),
    note_if_empty_2 = paste(
      "No pathway reached FDR < 0.05 in either discordant quadrant.",
      "The quadrants are small (18 and 7 proteins), so this is an absence of",
      "power rather than evidence that the discordant proteins share no biology."
    )
  ),
  ora_colors = ORA_QUAD_COLORS_F2,
  summary_quadrant_names = list(
    UU = "Concordant_Up", UU_slug = "concordant_up",
    DD = "Concordant_Down", DD_slug = "concordant_down",
    UD = "Discordant_YoungUp_OldDown", UD_slug = "discordant_y_up",
    DU = "Discordant_YoungDown_OldUp", DU_slug = "discordant_y_down"
  ),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  # Quadrant ORA needs all four corners populated; the discordant pair is
  # effectively empty here (UD = 0, DU = 1), so the supplementary block is off.
  supp = NULL
)
source("04_Figures/shared/comparison_panels/panel_E_rrho2.R")
panel_b <- pE_heat
n_shared_B <- n_shared

cfg <- c(engine_cfg, list(
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  panel_w = 102, panel_h = 140
))
source("04_Figures/shared/comparison_panels/panel_C_trajectory.R")
panel_c <- pC_trajectory

nes_scatter_cfg <- function(x, y, title, x_lab, y_lab, metric, interp,
                            slope, colors, draw_order, quadrants,
                            nudges = NULL, seed = 42) {
  list(
    fig_id = "F04", contrast_x = x, contrast_y = y, title = title,
    axis_x_label = x_lab, axis_y_label = y_lab,
    subtitle_metric = metric, subtitle_interpretation = interp,
    ref_slope = slope, panel_w = 120,
    rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
    sig_colors = colors,
    sig_draw_order = draw_order, quadrant_defs = quadrants,
    display_overrides = PATHWAY_LABELS, label_nudges = nudges,
    label_seed = seed
  )
}

# Empty until a render proves otherwise: box geometry drove the old nudges and
# the names are no longer boxed, so every one of them was re-derived from zero.
NUDGE_D <- tibble::tibble(
  pathway = character(), nudge_x = numeric(), nudge_y = numeric()
)

NUDGE_E <- NUDGE_D

cfg <- nes_scatter_cfg(
  "Training_Young", "Training_Old", "Pathway-Level Concordance (fGSEA)",
  "Training Young", "Training Old", "concordant",
  "positive ρ = shared pathway regulation across age groups", 1,
  SIG_COLORS_F2,
  c("Sig Old only", "Sig Young only", "Sig Both", "Interaction"),
  list(
    sig_both_label = "Sig Both", sig_x_label = "Sig Young only",
    sig_y_label = "Sig Old only",
    bg_red_1 = c(0, Inf, 0, Inf), bg_red_2 = c(-Inf, 0, -Inf, 0),
    bg_blue_1 = c(0, Inf, -Inf, 0), bg_blue_2 = c(-Inf, 0, 0, Inf),
    label_tr = "Concordant Up", color_tr = "#D6604D",
    label_tl = "Discordant", color_tl = "#4393C3",
    label_bl = "Concordant Down", color_bl = "#D6604D",
    label_br = "Discordant", color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q1 + q3
  ),
  nudges = NUDGE_D, seed = 7
)
source("04_Figures/shared/comparison_panels/panel_D_nes_scatter.R")
nes_legend <- cowplot::get_plot_component(pD, "guide-box-bottom", return_all = FALSE)
panel_d <- pD + theme(legend.position = "none")
rho_D <- as.numeric(nes_cor_all$estimate)
frac_D <- pw_conc_frac
file.rename(
  file.path(PNL_PNG, "MAIN_panel_D_nes_scatter.png"),
  file.path(PNL_PNG, "MAIN_panel_D_nes_concordance.png")
)
# The engine writes a fixed filename, so the second call would overwrite this.
file.rename(
  file.path(DAT, "panel_D", "nes_scatter.csv"),
  file.path(DAT, "panel_D", "nes_concordance.csv")
)

cfg <- nes_scatter_cfg(
  "Aging", "Training_Old", "Pathway-Level Reversal (fGSEA)",
  "Aging", "Training Old", "reversed",
  "negative ρ = training opposes aging effects", -1,
  SIG_COLORS_F3,
  c("Sig Training only", "Sig Aging only", "Sig Both"),
  list(
    sig_both_label = "Sig Both", sig_x_label = "Sig Aging only",
    sig_y_label = "Sig Training only",
    bg_blue_1 = c(0, Inf, -Inf, 0), bg_blue_2 = c(-Inf, 0, 0, Inf),
    bg_red_1 = c(0, Inf, 0, Inf), bg_red_2 = c(-Inf, 0, -Inf, 0),
    label_tr = "Exacerbated", color_tr = "#D6604D",
    label_tl = "Reversed", color_tl = "#4393C3",
    label_bl = "Exacerbated", color_bl = "#D6604D",
    label_br = "Reversed", color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q2 + q4
  ),
  nudges = NUDGE_E
)
source("04_Figures/shared/comparison_panels/panel_D_nes_scatter.R")
panel_e <- pD + theme(legend.position = "none")
rho_E <- as.numeric(nes_cor_all$estimate)
frac_E <- pw_conc_frac
file.rename(
  file.path(PNL_PNG, "MAIN_panel_D_nes_scatter.png"),
  file.path(PNL_PNG, "MAIN_panel_E_nes_reversal.png")
)
file.rename(
  file.path(DAT, "panel_D", "nes_scatter.csv"),
  file.path(DAT, "panel_D", "nes_reversal.csv")
)

cfg <- c(engine_cfg, list(
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  panel_w = 102, panel_h = 140
))
source("04_Figures/shared/comparison_panels/panel_fry_barcode.R")
panel_f <- pF_barcode

COMP_W <- 460
COMP_H <- 320
PRINT_SCALE <- COMP_W / 178
# The panels leave print_scale_apply.R's mutated FIG_* globals behind them, and
# composite_text_sizes multiplies by the canvas factor itself, so the constants
# have to be back at their style.R values or the titles scale twice. The panels
# are finished ggplot objects by now and keep the sizes they were built with.
source("04_Figures/shared/style.R")
# Authored at 460 mm and shrunk to the 178 mm print width, so the type scales
# up by the same factor to land at the shared 7 pt on the page.
f04_txt <- composite_text_sizes(COMP_W, 178)
TAG_SZ <- f04_txt$tag
TTL_SZ <- f04_txt$title
SUB_SZ <- f04_txt$subtitle

ttl <- c(
  A = "Training Concordance",
  B = "Threshold-Free Overlap",
  # C sits in a narrow 4-column box. At the 6 pt subtitle the title clips past
  # roughly 18 characters and the subtitle past roughly 25.
  C = "Blunted Response",
  D = "Pathway Concordance",
  E = "Pathway Reversal",
  F = "Signature Rank"
)
sub <- c(
  A = sprintf(
    # draw_label does not wrap and the eight-column box clips past about 80
    # characters. The fitted slopes live on the lines themselves instead.
    "N = %d, ρ %.2f | FDR %d: %.0f%% dir, ρ %.2f | Π %d: %.0f%%, ρ %.2f",
    n_total_A, r_spear,
    conc_fdr$n, 100 * conc_fdr$agree, conc_fdr$rho,
    conc_pi$n, 100 * conc_pi$agree, conc_pi$rho
  ),
  B = sweep_txt,
  C = "% young magnitude lost",
  D = sprintf("ρ = %.2f | %.0f%% concordant | no shared group mean", rho_D, frac_D * 100),
  E = sprintf("ρ = %.2f | %.0f%% reversed | pathway level only", rho_E, frac_E * 100),
  # draw_label does not wrap, and the six-column box runs out at about 54
  # characters, so the age wording is the short form here.
  F = "fry rotation test: younger sets in the older ranking"
)

layout <- paste(
  "##################",
  strrep("AAAAAAAABBBBBBCCCC\n", 6),
  "##################",
  strrep("DDDDDDEEEEEEFFFFFF\n", 6),
  "GGGGGGGGGGGGGGGGGG",
  sep = "\n"
)
layout <- paste(Filter(nzchar, strsplit(layout, "\n")[[1]]), collapse = "\n")

fig <- wrap_elements(full = panel_a) +
  # B's tile field sat 3.6 pt above panel A's border. The nudge is applied
  # here rather than in panel_E_rrho2.R because that engine also renders the
  # supplementary aging RRHO2, which must not move. B fills its cell height,
  # so a top margin lowers the top edge and shortens the square by the same
  # amount rather than translating it.
  wrap_elements(full = panel_b + theme(plot.margin = margin(1.27, 0, 0, 0, "mm"))) +
  wrap_elements(full = panel_c) +
  wrap_elements(full = panel_d) +
  wrap_elements(full = panel_e) +
  wrap_elements(full = panel_f) +
  wrap_elements(full = nes_legend) +
  plot_layout(
    # Panel C is given 18 pt more, taken from B's columns rather than A's. B is
    # a centred coord_fixed square in a cell 88 pt wider than the square needs,
    # so the width spent here was rendering as nothing and B stays
    # height-limited at the same size. Taking it from A instead would reflow
    # its ggrepel labels and push its bar labels further past their bars.
    design = layout, widths = c(rep(1, 8), rep(0.958, 6), rep(1.063, 4)),
    # Bands 1 and 8 carry the title and subtitle for the row below them. They
    # are fractions of the total, so trimming the middle band shrinks the top
    # one too, and the top band had to grow from 7 to hold its ground. The
    # middle band was 12, which left 26 pt of dead canvas between the D/E/F
    # subtitles and the plots they name; at 9 that gap is about 5 pt.
    heights = c(8, rep(13, 6), 9, rep(13, 6), 6)
  )

# Hand-tuned to the column boundaries, with no link back to the layout, so a
# change to widths has to be carried here by hand. C's cell left edge moved
# 18 pt left with the widths change above; B's did not move.
X_TAG <- c(A = 0.004, B = 0.448, C = 0.768, D = 0.004, E = 0.334, F = 0.659)
# The top row's tags started 0.6 pt from the page edge, against the bottom
# row's 2.1 pt below its band. Lowered to match that relationship and to put
# a readable gap between the subtitles and the panel tops at 52.4.
Y_TAG <- c(A = 0.987, B = 0.987, C = 0.987, D = 0.518, E = 0.518, F = 0.518)
X_TTL <- 0.030
SUB_OFFSET <- 0.019

composite_final <- Reduce(
  function(p, k) {
    p +
      draw_label(k, x = X_TAG[[k]], y = Y_TAG[[k]], size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
      draw_label(ttl[[k]],
        x = X_TAG[[k]] + X_TTL, y = Y_TAG[[k]], size = TTL_SZ,
        fontface = "bold", hjust = 0, vjust = 1
      ) +
      draw_label(sub[[k]],
        x = X_TAG[[k]] + X_TTL, y = Y_TAG[[k]] - SUB_OFFSET, size = SUB_SZ,
        fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40"
      )
  },
  names(ttl),
  init = ggdraw(fig)
)

ggsave(file.path(OUT_PDF, "MAIN_F04_composite.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
ggsave(file.path(OUT_PNG, "MAIN_F04_composite.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)

message("F04 composite (6-panel, A/B/C over D/E/F) saved")
