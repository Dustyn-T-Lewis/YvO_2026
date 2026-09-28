#!/usr/bin/env Rscript
# S5 Table, sheet SUPP_threshold_sens: diagnostic for main panel B, showing the
# concordance pattern is stable across significance thresholds (Pi < 0.05, FDR < 0.05, FDR < 0.10, nominal p < 0.05).

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

DAT <- "04_Figures/F04/c_data/panel_supp"
RPT <- "04_Figures/F04/b_reports/supp/panels"
for (d in c(DAT, RPT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)

base_df <- dep_df |>
  transmute(gene,
            logFC_TY   = logFC_Training_Young,
            logFC_TO   = logFC_Training_Old,
            pi_TY      = pi_score_Training_Young,
            pi_TO      = pi_score_Training_Old,
            fdr_TY     = adj.P.Val_Training_Young,
            fdr_TO     = adj.P.Val_Training_Old,
            nom_TY     = P.Value_Training_Young,
            nom_TO     = P.Value_Training_Old) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO))

assign_quad <- function(lfc_ty, lfc_to) {
  case_when(
    lfc_ty > 0 & lfc_to > 0 ~ "Concordant Up",
    lfc_ty < 0 & lfc_to < 0 ~ "Concordant Down",
    TRUE                     ~ "Discordant")
}

thresholds <- list(
  "Π < 0.05"      = function(d) d |> filter(pi_TY < 0.05 | pi_TO < 0.05),
  "FDR < 0.05"     = function(d) d |> filter(fdr_TY < 0.05 | fdr_TO < 0.05),
  "FDR < 0.10"     = function(d) d |> filter(fdr_TY < 0.10 | fdr_TO < 0.10),
  "Nom. p < 0.05" = function(d) d |> filter(nom_TY < 0.05 | nom_TO < 0.05)
)

sens_df <- bind_rows(lapply(names(thresholds), \(thr_name) {
  sig_sub <- thresholds[[thr_name]](base_df) |>
    mutate(quadrant = assign_quad(logFC_TY, logFC_TO))
  counts <- sig_sub |> count(quadrant, name = "n")
  counts$threshold <- thr_name
  counts$n_total <- nrow(sig_sub)
  counts
})) |>
  mutate(threshold = factor(threshold, levels = names(thresholds)),
         quadrant  = factor(quadrant, levels = c("Concordant Up",
                                                  "Concordant Down",
                                                  "Discordant")))

write_csv(sens_df, file.path(DAT, "SUPP_threshold_sensitivity.csv"))
message("Threshold sensitivity:\n", paste(capture.output(print(sens_df)), collapse = "\n"))

QUAD_COLS <- c("Concordant Up" = "#E57373",
               "Concordant Down" = "#64B5F6",
               "Discordant" = "#FFB74D")

pS_thresh <- ggplot(sens_df, aes(x = threshold, y = n, fill = quadrant)) +
  geom_col(position = "stack", width = 0.65,
           color = "white", linewidth = 0.3) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            size = BASE_COUNT, fontface = "bold", color = "white") +
  scale_fill_manual(values = QUAD_COLS, name = "Quadrant") +
  labs(title = "Concordance Pattern: Threshold Sensitivity",
       subtitle = "Protein counts per quadrant across significance criteria",
       x = NULL, y = "Significant proteins") +
  FIG_THEME

PW <- 89; PH <- 70
ggsave(file.path(RPT, "S5T_threshold_sens.png"), pS_thresh,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT, "S5T_threshold_sens.pdf"), pS_thresh,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("SUPP Panel C (threshold sensitivity) done")

invisible(pS_thresh)
