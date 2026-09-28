#!/usr/bin/env Rscript
# S5a Figure B: diagnostic for main panel B, the GO Slim category breakdown by
# concordance quadrant.

setwd(here::here())

pacman::p_load(dplyr, tidyr, tibble, stringr, readr, ggplot2, patchwork, cowplot)

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

source("04_Figures/shared/go_slim_categories.R")
pacman::p_load(readxl)

DAT <- "04_Figures/F04/c_data/panel_supp"
RPT <- "04_Figures/F04/b_reports/supp/panels"
for (d in c(DAT, RPT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pattern_df <- read_csv(DEP_RESULTS, show_col_types = FALSE) |>
  filter(!is.na(logFC_Training_Young), !is.na(logFC_Training_Old)) |>
  filter(
    pi_score_Training_Young < 0.05 | pi_score_Training_Old < 0.05 |
      pi_score_Interaction < 0.05
  ) |>
  mutate(quadrant = case_when(
    logFC_Training_Young > 0 & logFC_Training_Old > 0 ~ "Concordant Up",
    logFC_Training_Young < 0 & logFC_Training_Old < 0 ~ "Concordant Down",
    TRUE ~ "Discordant"
  )) |>
  as.data.frame()

slim_merged <- assign_go_slim_consolidated(fg_genes = pattern_df$gene, all_genes = pattern_df$gene) |>
  left_join(transmute(pattern_df, gene, quadrant), by = "gene") |>
  filter(!is.na(quadrant), !is.na(consolidated))

cat_quad <- slim_merged |>
  count(consolidated, quadrant, name = "n") |>
  mutate(quadrant = factor(quadrant,
                           levels = c("Concordant Up", "Concordant Down",
                                      "Discordant")))

cat_order <- cat_quad |>
  group_by(consolidated) |>
  summarise(total = sum(n), .groups = "drop") |>
  arrange(total) |>
  pull(consolidated)
cat_order <- c("Other", setdiff(cat_order, "Other"))

cat_quad <- cat_quad |>
  mutate(consolidated = factor(consolidated, levels = cat_order))

write_csv(cat_quad, file.path(DAT, "SUPP_goslim_distribution.csv"))
message("GO Slim distribution:\n", paste(capture.output(print(cat_quad)), collapse = "\n"))

QUAD_COLS <- c("Concordant Up" = "#E57373",
               "Concordant Down" = "#64B5F6",
               "Discordant" = "#FFB74D")

pS_goslim <- ggplot(cat_quad, aes(x = n, y = consolidated, fill = quadrant)) +
  geom_col(position = "stack", width = 0.7,
           color = "white", linewidth = 0.3) +
  scale_fill_manual(values = QUAD_COLS, name = "Quadrant") +
  labs(title = "GO Slim Category Distribution",
       subtitle = "Protein counts per GO Slim category by concordance quadrant",
       x = "Protein count", y = NULL) +
  FIG_THEME

PW <- 89; PH <- 70
ggsave(file.path(RPT, "S5a_B_goslim_bars.png"), pS_goslim,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT, "S5a_B_goslim_bars.pdf"), pS_goslim,
       width = PW, height = PH, units = "mm", device = get_pdf_device())

message("SUPP Panel D (GO Slim bars) done")

invisible(pS_goslim)
