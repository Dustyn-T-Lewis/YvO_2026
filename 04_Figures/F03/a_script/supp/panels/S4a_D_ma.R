#!/usr/bin/env Rscript
# S4a Figure D: MA plots, all four contrasts, coloured by Pi direction.

setwd(here::here())
source("04_Figures/F03/a_script/supp/panels/_distribution.R", local = TRUE)

RPT <- "04_Figures/F03/b_reports/supp/panels"
DAT <- "04_Figures/F03/c_data/supp"
for (d in c(RPT, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
per_contrast <- read_contrasts()

ma_df <- bind_rows(lapply(CTRS, \(ctr) {
  per_contrast[[ctr]] |>
    transmute(contrast = ctr, gene, average_intensity, logFC,
              sig_pi, direction = case_when(
                sig_pi ==  1 ~ "Up", sig_pi == -1 ~ "Down", TRUE ~ "NS"))
})) |>
  mutate(contrast  = factor(contrast, levels = CTRS),
         direction = factor(direction, levels = c("Up", "Down", "NS")))

n_dep <- ma_df |> filter(direction != "NS") |> count(contrast, name = "n_dep")

p <- ggplot(ma_df, aes(average_intensity, logFC, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(data = \(d) filter(d, direction == "NS"), alpha = 0.25, size = 0.6) +
  geom_point(data = \(d) filter(d, direction != "NS"), alpha = 0.85, size = 0.9) +
  geom_text(data = n_dep, aes(x = Inf, y = Inf, label = sprintf("Π: %d", n_dep)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 2.2,
            fontface = "bold", color = "grey20") +
  scale_color_manual(values = DIR_COLORS, name = "Π < 0.05") +
  facet_wrap(~contrast, ncol = 2, labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "MA plots", x = "Mean log2 intensity", y = "logFC", tag = "d") +
  FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank(),
                    legend.position = "top", legend.key.size = unit(3, "mm"))
write_csv(ma_df, file.path(DAT, "panel_C_ma.csv"))
ggsave(file.path(RPT, "S4a_D_ma.png"), p, width = 89, height = 75,
       units = "mm", dpi = 300)
ggsave(file.path(RPT, "S4a_D_ma.pdf"), p, width = 89, height = 75,
       units = "mm", device = get_pdf_device())

invisible(p)
