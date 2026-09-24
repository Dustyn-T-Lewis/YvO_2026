# Pathway NES scatter config, shared by panels D and E: the label shortenings,
# the engine config builder and the label nudges.

PNL_PNG <- "04_Figures/F04/b_reports/panels"
PNL_PDF <- PNL_PNG
DAT <- "04_Figures/F04/c_data"

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

# The engine writes fixed file names, and D and E both run it, so each panel
# moves its outputs to its own names before the other can overwrite them.
rename_nes_outputs <- function(name, csv) {
  for (ext in c("png", "pdf")) {
    file.rename(
      file.path(PNL_PNG, paste0("MAIN_panel_D_nes_scatter.", ext)),
      file.path(PNL_PNG, paste0(name, ".", ext))
    )
  }
  file.rename(
    file.path(PNL_PNG, "MAIN_panel_D_legend.png"),
    file.path(PNL_PNG, paste0(substr(name, 1, 1), "_legend.png"))
  )
  file.rename(
    file.path(DAT, "panel_D", "nes_scatter.csv"),
    file.path(DAT, "panel_D", csv)
  )
}
