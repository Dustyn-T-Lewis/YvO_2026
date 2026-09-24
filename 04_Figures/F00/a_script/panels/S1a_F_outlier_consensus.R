#!/usr/bin/env Rscript
# S1a Figure F: outlier detection consensus.

setwd(here::here())
source("04_Figures/F00/a_script/panels/_inputs.R", local = TRUE)

od <- int_norm$outlier_diag |>
  mutate(
    age = factor(ifelse(grepl("^O", prefix), "Old", "Young"),
      levels = c("Young", "Old")
    ),
    status = case_when(
      consensus_outlier ~ "Outlier",
      n_flags > 0 ~ "Flagged",
      TRUE ~ "Clean"
    )
  )

# geom_jitter draws from the RNG, so without this the panel's point positions
# move on every run and the figure is not reproducible.
set.seed(42)
pF <- ggplot(od, aes(mahal_dist, n_flags, color = status, shape = Timepoint)) +
  geom_jitter(width = 0.04, height = 0.12, size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c(
    Clean = "grey55", Flagged = "#E6A100",
    Outlier = "#B2182B"
  ), name = NULL) +
  scale_shape_manual(values = SHAPE_TP, guide = "none") +
  scale_y_continuous(breaks = 0:4) +
  labs(
    x = "Mahalanobis distance", y = "QC flags", tag = "F",
    title = "Outlier detection consensus",
    subtitle = sprintf(
      "%d/%d removed: %s",
      sum(od$consensus_outlier), nrow(od),
      paste(int_norm$outlier_ids, collapse = ", ")
    )
  ) +
  FIG_THEME +
  theme(legend.position = "top", legend.key.size = unit(2.5, "mm"))

PNL_F <- file.path(PNL, "S1a_F_outlier_consensus")
ggsave(paste0(PNL_F, ".png"), pF,
  width = PW, height = PH, units = "mm", dpi = 300
)
# The jitter seed is drawn from the RNG each time the plot is rendered. Seeding
# again gives the PDF the PNG's points and leaves the RNG where the PNG left
# it, which is where S1a.R's composite renders have always started from.
set.seed(42)
ggsave(paste0(PNL_F, ".pdf"), pF,
  width = PW, height = PH, units = "mm", device = get_pdf_device()
)
saveRDS(od, file.path(SHEETS, "panel_F.rds"))

invisible(pF)
