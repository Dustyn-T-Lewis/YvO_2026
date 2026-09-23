# Shared helpers for 90_stitch_figure.R scripts: safe CSV reads, workbook assembly, cleanup.

pacman::p_load(openxlsx, readr, readxl)

add_sheet <- function(wb, name, data) {
  # A zero-row frame writes a sheet with no header at all, which reads as a
  # broken export rather than an empty result.
  if (nrow(data) == 0) {
    data <- data.frame(
      note = "No rows: the analysis backing this sheet returned no results.",
      stringsAsFactors = FALSE
    )
  }
  addWorksheet(wb, name)
  writeData(wb, name, data)
  hs <- createStyle(textDecoration = "bold")
  addStyle(wb, name, hs, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
  cat(sprintf("    + %s: %d x %d\n", name, nrow(data), ncol(data)))
}

safe_read <- function(path) {
  if (file.exists(path)) {
    as.data.frame(read_csv(path))
  } else {
    cat(sprintf("    SKIP (not found): %s\n", path))
    NULL
  }
}

build_workbook <- function(out_file, title = NULL, description = NULL,
                           overview_df = NULL, sheet_specs) {
  # Overview first, always. These workbooks are the journal's supplementary
  # tables; the largest carries 22 sheets and F00's are named panel_A..panel_N,
  # so without an index the only map is the R script that wrote them. The
  # overview_df the callers build used to be accepted and discarded.
  #
  # Sheets that safe_read() skipped are left out, so the index cannot promise a
  # sheet the workbook does not contain.
  wb <- createWorkbook()
  addWorksheet(wb, "Overview")

  written <- list()
  for (spec in sheet_specs) {
    df <- if (!is.null(spec$df)) spec$df else safe_read(spec$path)
    if (!is.null(df)) {
      add_sheet(wb, spec$name, df)
      written[[length(written) + 1L]] <- data.frame(
        Sheet = spec$name, Rows = nrow(df), Columns = ncol(df)
      )
    }
  }
  index <- do.call(rbind, written)

  described <- if (!is.null(overview_df) && nrow(overview_df)) {
    overview_df$Description[match(index$Sheet, overview_df$Sheet)]
  } else {
    rep(NA_character_, nrow(index))
  }
  index$Description <- ifelse(is.na(described), "", described)

  head_rows <- c(title, description)
  if (length(head_rows)) {
    writeData(wb, "Overview", data.frame(head_rows), startRow = 1, colNames = FALSE)
    addStyle(wb, "Overview", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  }
  start <- length(head_rows) + if (length(head_rows)) 2L else 1L
  writeData(wb, "Overview", index, startRow = start)
  addStyle(wb, "Overview", createStyle(textDecoration = "bold"),
    rows = start, cols = seq_len(ncol(index)), gridExpand = TRUE
  )
  setColWidths(wb, "Overview", cols = seq_len(ncol(index)), widths = "auto")

  saveWorkbook(wb, out_file, overwrite = TRUE)
  cat(sprintf(
    "  Saved: %s (%.0f KB, %d sheets + Overview)\n\n",
    out_file, file.size(out_file) / 1e3, nrow(index)
  ))
}

# Cross-figure readers: prefer these over reading raw CSVs from another figure's c_data/.
read_sheet_df <- function(xlsx, sheet) {
  stopifnot("supplementary workbook missing" = file.exists(xlsx))
  as.data.frame(readxl::read_excel(xlsx, sheet = sheet))
}

read_matrix_sheet <- function(xlsx, sheet, row_col = "sample_id") {
  df <- read_sheet_df(xlsx, sheet)
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df[[row_col]]
  mat
}

read_vector_sheet <- function(xlsx, sheet) read_sheet_df(xlsx, sheet)[[1]]

matrix_to_df <- function(mat, row_col = "sample_id") {
  df <- data.frame(
    rn = rownames(mat),
    as.data.frame(mat, check.names = FALSE),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  names(df)[1] <- row_col
  df
}

# Remove per-panel CSV intermediates after the workbook is written.
# preserve_patterns: regexes for paths that must not be deleted (upstream stages, shared cache).
cleanup_after_workbook <- function(sheet_specs,
                                   extra_subdirs = character(),
                                   extra_files = character(),
                                   preserve_patterns = UPSTREAM_PREFIXES) {
  is_preserved <- function(path) {
    rel <- sub(paste0("^", "", "/?"), "", path)
    any(vapply(preserve_patterns, function(p) grepl(p, path) || grepl(p, rel), logical(1)))
  }
  removed <- 0L
  preserved <- 0L
  for (sp in sheet_specs) {
    if (!is.null(sp$path) && file.exists(sp$path)) {
      if (is_preserved(sp$path)) {
        preserved <- preserved + 1L
      } else {
        unlink(sp$path)
        removed <- removed + 1L
      }
    }
  }
  for (sd in extra_subdirs) {
    if (dir.exists(sd)) {
      unlink(sd, recursive = TRUE)
      removed <- removed + 1L
    }
  }
  for (f in extra_files) {
    if (file.exists(f) && !is_preserved(f)) {
      unlink(f)
      removed <- removed + 1L
    }
  }
  cat(sprintf(
    "  cleanup: removed %d intermediate(s); preserved %d upstream/shared path(s)\n",
    removed, preserved
  ))
}
