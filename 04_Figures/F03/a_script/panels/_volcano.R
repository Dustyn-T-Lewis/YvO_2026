# Volcano ring for one contrast, shared by panels A-D. Saves the standalone
# panel and the ring's term table, and returns the plot with its tag.

pacman::p_load(readr, dplyr, ggplot2)

source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")
source("04_Figures/shared/build_fgsea_cache.R")

volcano_panel <- function(contrast, title, subtitle, tag, name) {
  VW <- 89
  VH <- 89 # half of 178mm double-col -> no journal scaling
  RPT <- "04_Figures/F03/b_reports/panels"
  DAT <- file.path("04_Figures/F03/c_data", paste0("panel_", tag))
  dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

  dep_df <- read_csv("03_DEP/c_data/03_combined_results.csv", show_col_types = FALSE)
  fgsea_all <- read_csv("04_Figures/shared/fgsea_tstat_all_v2.csv", show_col_types = FALSE)

  pi_col <- paste0("pi_score_", contrast)
  n_dep <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  n_path <- sum(!is.na(fgsea_all$padj) & fgsea_all$padj < 0.05 &
    fgsea_all$contrast == contrast &
    fgsea_all$database %in% c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome"))
  enriched_sub <- sprintf("%s | %d DEPs, %d pathways", subtitle, n_dep, n_path)

  top_terms <- select_ring_terms(fgsea_all, contrast)
  ring_data <- build_ring_with_gaps(top_terms, contrast, fgsea_all)
  max_arc <- max(ring_data$arc_r1_var, na.rm = TRUE)
  adaptive_gap <- 0.7 + 0.3 * (max_arc - 4.8) / 1.6

  p <- make_volcano_ring(
    de_df              = dep_df,
    go_df              = fgsea_all,
    contrast           = contrast,
    title              = NULL,
    contrast_title     = title,
    contrast_subtitle  = enriched_sub,
    ring_data_override = ring_data,
    label_size         = 2.35,
    label_gap          = adaptive_gap,
    title_size         = 5,
    subtitle_size      = 3.5,
    point_size         = 0.5,
    point_alpha        = 0.55,
    count_label_size   = scale_text(BASE_COUNT, VW) + 0.4,
    count_y_mult       = 0.75,
    count_x_mult       = 0.85,
    bg_color           = unname(CONTRAST_COLORS[contrast]),
    bg_alpha           = 0.20,
    show_legend        = FALSE
  ) + labs(tag = tag)

  ggsave(file.path(RPT, paste0(name, ".png")), p,
    width = VW, height = VH, units = "mm", dpi = 300
  )
  ggsave(file.path(RPT, paste0(name, ".pdf")), p,
    width = VW, height = VH, units = "mm", device = get_pdf_device()
  )

  ring_out <- attr(p, "ring_data")
  if (!is.null(ring_out) && nrow(ring_out) > 0) {
    dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
    write_csv(
      ring_out |> dplyr::select(
        -gene_list,
        -any_of(c(
          "term_idx", "start_deg", "end_deg", "mid_deg",
          "start_rad", "end_rad", "mid_rad", "arc_r1_var"
        ))
      ),
      file.path(DAT, "ring_terms.csv")
    )
  }

  message(sprintf("F03 panel %s (%s) done", tag, contrast))
  invisible(p)
}
