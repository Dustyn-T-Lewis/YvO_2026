#!/usr/bin/env Rscript
# F04 — Training Concordance: Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main panels (composite + xlsx + cleanup).

setwd(here::here())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT <- "04_Figures/F04/c_data"

# Supp panels first (CSVs needed for xlsx)
source("04_Figures/F04/a_script/02_supp_panels.R")

# Main panels + composite
source("04_Figures/F04/a_script/01_main_panels.R")

# Build supplementary xlsx
message("F04 supplementary workbook")
enrichment_blunting_df <- read.csv(file.path(DAT, "panel_supp", "enrichment_blunting.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
f04_specs <- list(
  list(name = "panel_A_ora_quadrant", path = file.path(DAT, "panel_A", "ora_quadrant.csv")),
  list(name = "panel_B_rrho2_summary", path = file.path(DAT, "panel_E", "rrho2_summary.csv")),
  list(name = "panel_B_rrho2_hotspot", path = file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "panel_B_rrho2_ora_concord", path = file.path(DAT, "panel_E", "rrho2_ora_concordant.csv")),
  list(name = "panel_B_rrho2_ora_discord", path = file.path(DAT, "panel_E", "rrho2_ora_discordant.csv")),
  list(name = "panel_C_trajectory", path = file.path(DAT, "panel_C_trajectory.csv")),
  list(name = "panel_C_trajectory_prot", path = file.path(DAT, "panel_C_trajectory_proteins.csv")),
  list(name = "panel_D_nes_concordance", path = file.path(DAT, "panel_D", "nes_concordance.csv")),
  list(name = "panel_E_nes_reversal", path = file.path(DAT, "panel_D", "nes_reversal.csv")),
  list(name = "panel_F_fry_barcode", path = file.path(DAT, "panel_F_fry_barcode.csv")),
  list(name = "SUPP_enrichment_blunting", df = enrichment_blunting_df),
  list(name = "SUPP_ora_dedup", path = file.path(DAT, "panel_supp", "SUPP_ora_dedup_sensitivity.csv")),
  list(name = "SUPP_rho_bootstrap", path = file.path(DAT, "panel_supp", "SUPP_rho_bootstrap.csv")),
  list(name = "SUPP_threshold_sens", path = file.path(DAT, "panel_supp", "SUPP_threshold_sensitivity.csv")),
  list(name = "SUPP_goslim_bars", path = file.path(DAT, "panel_supp", "SUPP_goslim_distribution.csv")),
  list(name = "SUPP_fry_leading", path = file.path(DAT, "panel_supp", "SUPP_fry_leading_edge.csv")),
  list(name = "SUPP_rrho2_aging_summary", path = file.path(DAT, "panel_supp", "panel_E", "rrho2_summary.csv")),
  list(name = "SUPP_rrho2_aging_hotspot", path = file.path(DAT, "panel_supp", "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "SUPP_young_dep_heatmap", path = file.path(DAT, "panel_supp", "SUPP_young_dep_heatmap.csv")),
  list(name = "SUPP_young_dep_ora_conc", path = file.path(DAT, "panel_supp", "SUPP_young_dep_ora_concordant.csv")),
  list(name = "SUPP_young_dep_discord", path = file.path(DAT, "panel_supp", "SUPP_young_dep_discordant.csv")),
  list(name = "SUPP_young_dep_ora", path = file.path(DAT, "panel_supp", "SUPP_young_dep_ora.csv")),
  list(name = "SUPP_concordance_magnitude", path = file.path(DAT, "panel_supp", "SUPP_concordance_magnitude.csv")),
  list(name = "SUPP_cat_depth", path = file.path(DAT, "panel_supp", "SUPP_cat_depth.csv")),
  list(name = "SUPP_coupling_null", path = file.path(DAT, "panel_supp", "SUPP_coupling_null.csv"))
)
build_workbook(
  file.path(DAT, "F04_supplementary.xlsx"),
  title = "S5 Table \u2014 training-response concordance",
  description = "Training-concordance diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_rrho2_summary", "panel_B_rrho2_hotspot",
      "panel_B_rrho2_ora_concord", "panel_B_rrho2_ora_discord",
      "panel_C_trajectory", "panel_C_trajectory_prot",
      "panel_D_nes_concordance",
      "panel_E_nes_reversal",
      "panel_F_fry_barcode",
      "SUPP_enrichment_blunting",
      "SUPP_ora_dedup", "SUPP_rho_bootstrap", "SUPP_threshold_sens",
      "SUPP_goslim_bars", "SUPP_fry_leading",
      "SUPP_rrho2_aging_summary", "SUPP_rrho2_aging_hotspot",
      "SUPP_concordance_magnitude", "SUPP_cat_depth", "SUPP_coupling_null",
      "SUPP_young_dep_heatmap", "SUPP_young_dep_ora_conc",
      "SUPP_young_dep_discord", "SUPP_young_dep_ora"
    ),
    Description = c(
      "Panel A: ORA by training-concordance scatter quadrant",
      "Panel B: RRHO2 quadrant summary (max -log10p per quadrant)",
      "Panel B: RRHO2 hotspot genes per quadrant",
      "Panel B: ORA on RRHO2 concordant quadrant genes",
      "Panel B: ORA on RRHO2 discordant quadrant genes",
      "Panel C: response-magnitude compression summary (concordance, sign changes, proteome-wide ratio)",
      "Panel C: the young-responsive proteins drawn in the panel, with both training log2FCs",
      "Panel D: NES scatter (Training_Young vs Training_Old) per pathway, Spearman + Fisher Z",
      "Panel E: NES scatter (Aging vs Training_Old) per pathway, Spearman + Fisher Z",
      "Panel F: fry rotation test for Training_Young signatures against the Training_Old ranking",
      "SUPP: blunting enrichment (Pi-score-weighted ORA on Training_Old DEPs)",
      "SUPP: ORA dedup sensitivity across Jaccard cutoffs",
      "SUPP: Spearman rho bootstrap (1000 reps, 95% CI)",
      "SUPP: Protein counts per quadrant across significance thresholds",
      "SUPP: GO Slim category distribution by concordance quadrant",
      "SUPP: Top 20 fry driving proteins by |t-stat| in Training Old",
      "SUPP: Aging x Training_Old RRHO2 quadrant summary, descriptive (contrasts share Old_Pre; top-k columns are diagnostics, not a test)",
      "SUPP: Aging x Training_Old RRHO2 hotspot genes per quadrant",
      "SUPP panel C: sign agreement across quintiles of min(|log2FC|)",
      "SUPP: corner enrichment against list depth (CAT curve). The diagnostic behind main panel B, whose corner test reads one point off this curve at the top 10%",
      "SUPP panel D: observed Aging-vs-Training_Old r against a coupling-preserving null",
      "S5b Figure: per-participant training change for the 135 younger-adult FDR proteins",
      "S5b Figure: over-representation for the 110 proteins that move the same way in both age groups",
      "S5b Figure: the 25 proteins that move the opposite way, with both training log2FCs",
      "S5b Figure: over-representation for the increased and the decreased proteins separately"
    ),
    stringsAsFactors = FALSE
  ),
  sheet_specs = f04_specs
)
cleanup_after_workbook(f04_specs,
  extra_subdirs = c(
    file.path(DAT, "panel_A"),
    file.path(DAT, "panel_D"),
    file.path(DAT, "panel_E"),
    file.path(DAT, "panel_supp")
  )
)

remaining <- list.files(DAT, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("F04 complete")
