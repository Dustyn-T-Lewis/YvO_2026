#!/usr/bin/env Rscript
# S5 Table, sheet SUPP_ora_dedup: ORA Dedup Sensitivity
# Diagnostic for main Panel A: ORA results do not depend on the Jaccard dedup cutoff.

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

pdf_device <- get_pdf_device()

BASE <- "04_Figures/F04"
DAT  <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv",
                   show_col_types = FALSE)

scatter_df <- dep_df |>
  transmute(gene,
            logFC_TY = logFC_Training_Young,
            logFC_TO = logFC_Training_Old) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO)) |>
  mutate(quadrant = case_when(
    logFC_TY > 0 & logFC_TO > 0 ~ "Concordant Up",
    logFC_TY < 0 & logFC_TO < 0 ~ "Concordant Down",
    logFC_TY > 0 & logFC_TO < 0 ~ "Discordant (Y Up / O Down)",
    TRUE                         ~ "Discordant (Y Down / O Up)"))

universe <- scatter_df$gene

# Pathway collection (same as panel_A_ORA.R)
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# Sweep the EnrichmentMap combined coefficient, the rule the analysis uses.
# 1.0 collapses only identical sets, so it reads as the no-collapsing reference.
cutoffs <- c(0.25, 0.375, 0.5, 1.0)
target_quads <- c("Concordant Up", "Concordant Down")

results <- list()
for (quad in target_quads) {
  genes <- scatter_df$gene[scatter_df$quadrant == quad]
  if (length(genes) < 5) next

  for (jc in cutoffs) {
    message(sprintf("  ORA: %s | EM coefficient = %.3f", quad, jc))
    ora_res <- tryCatch({
      run_ora_deduplicated(
        genes = genes, universe = universe, pathways = pw_collection,
        em_cutoff = jc, min_size = 15, max_size = 500, padj_cutoff = 0.05
      )
    }, error = function(e) { message("    error: ", e$message); tibble() })

    n_sig <- if (nrow(ora_res) > 0) {
      sum(ora_res$dedup_status == "kept", na.rm = TRUE)
    } else {
      0L
    }
    results[[length(results) + 1]] <- tibble(
      quadrant = quad, em_cutoff = jc,
      n_representatives = n_sig, n_genes = length(genes))
  }
}

sens_df <- bind_rows(results)
write_csv(sens_df, file.path(DAT, "SUPP_ora_dedup_sensitivity.csv"))
message("Sensitivity data:\n", paste(capture.output(print(sens_df)), collapse = "\n"))

sens_df <- sens_df |>
  mutate(cutoff_label = factor(sprintf("%.3g", em_cutoff),
                               levels = sprintf("%.3g", cutoffs)))

pS_ora_dedup <- ggplot(sens_df, aes(x = quadrant, y = n_representatives,
                                     fill = cutoff_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n_representatives),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = BASE_COUNT, fontface = "bold") +
  scale_fill_manual(values = c("0.25" = "grey80", "0.375" = "grey55",
                               "0.5" = "grey35", "1" = "grey15"),
                    name = "EM coefficient") +
  labs(title = "ORA Redundancy Sensitivity",
       subtitle = "Representative pathways (FDR < 0.05) by EnrichmentMap cutoff",
       x = NULL, y = "Representative pathways") +
  FIG_THEME +
  theme(legend.key.size = unit(1.8, "mm"))

RPT_PNG <- file.path(BASE, "b_reports", "supp", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "S5T_ora_dedup.png"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "S5T_ora_dedup.pdf"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("SUPP Panel A (ORA dedup) done")

invisible(pS_ora_dedup)
