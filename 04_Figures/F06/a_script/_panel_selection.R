#!/usr/bin/env Rscript
# Which (comparison x module) cells the main figure draws.
#
# Four scripts need this list: the two ROC panels and the two LOSO sensitivity
# analyses, whose S7 Table sheets are meant to cover what the figure shows.
# Each used to rebuild it from its own copy of the same arrange-and-slice, so a
# change to the figure left the tables describing a different set of cells.
#
# The figure is split by comparison rather than ranked across all of them: one
# panel asks whether a module separates the age groups, the other whether it
# separates pre from post. A cell earns its place by clearing p < 0.05 in its
# own comparison, so neither panel carries a cell that answers nothing.

# The training panel is built from the modules that respond in the younger
# group, each shown again in the older group. Nothing reaches p < 0.05 in the
# older group, so selecting on it would empty the panel; the younger fit is
# what the older cell is being compared against.
main_panel_cells <- function(summ) {
  age <- summ |>
    dplyr::filter(row == "Age", perm_p < 0.05) |>
    dplyr::arrange(dplyr::desc(auc)) |>
    dplyr::mutate(panel = "age", stratum = NA_character_)

  keep <- summ |>
    dplyr::filter(row == "Pre vs Post (Young)", perm_p < 0.05) |>
    dplyr::arrange(dplyr::desc(auc)) |>
    dplyr::pull(module)

  training <- summ |>
    dplyr::filter(
      row %in% c("Pre vs Post (Young)", "Pre vs Post (Old)"),
      module %in% keep
    ) |>
    dplyr::mutate(
      panel = "training",
      stratum = ifelse(grepl("Young", row), "Younger", "Older"),
      module = factor(module, levels = keep)
    ) |>
    dplyr::arrange(dplyr::desc(stratum), module) |>
    dplyr::mutate(module = as.character(module))

  dplyr::bind_rows(age, training) |>
    dplyr::mutate(
      sig = dplyr::case_when(
        !is.na(q_bh) & q_bh < 0.05 ~ "q<.05",
        perm_p < 0.05 ~ "p<.05",
        TRUE ~ "ns"
      ),
      # Thinner than the supplement's 2.4 pt: these cells are a third the size.
      border_color = ifelse(sig == "ns", "grey80", "black"),
      border_lty = ifelse(sig == "p<.05", "dashed", "solid"),
      border_lw = ifelse(sig == "ns", 0.2, 0.8),
      cell_label = ifelse(
        panel == "age",
        stringr::str_to_title(module),
        paste0(stringr::str_to_title(module), " · ", stratum)
      )
    )
}
