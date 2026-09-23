# An Overview sheet for the three stage workbooks.
#
# The seven figure workbooks get theirs from build_workbook() in
# figure_supplement_helpers.R. The stage workbooks are assembled sheet by sheet
# instead, and shipped as S8, S9 and S10 Table with no index at all: a reviewer
# opening one met `normalized_matrix` and `mnar_audit` with nothing to read them
# against.
#
# Sheet names stay as they are. They are keys the figure scripts read by name,
# so the index explains them rather than renaming them.

pacman::p_load(openxlsx)

add_overview <- function(wb, title, description, entries) {
  # A description for a sheet that was not written is dropped: one sheet is
  # conditional on its analysis returning rows, and an index must not promise a
  # sheet the workbook does not contain. A sheet with no description is an
  # error, because that is the failure this exists to prevent.

  # The index is rebuilt by whichever script writes the workbook last, so a
  # script that appends to an already-indexed workbook has to be able to run.
  if ("Overview" %in% names(wb)) removeWorksheet(wb, "Overview")

  present <- names(wb)
  undescribed <- setdiff(present, names(entries))
  if (length(undescribed)) {
    stop("add_overview: sheet with no description: ",
         paste(undescribed, collapse = ", "))
  }

  index <- data.frame(
    Sheet = present,
    Rows = vapply(present, \(s) nrow(readWorkbook(wb, s)), integer(1)),
    Columns = vapply(present, \(s) ncol(readWorkbook(wb, s)), integer(1)),
    Contents = unname(entries[present]),
    stringsAsFactors = FALSE
  )

  addWorksheet(wb, "Overview")
  writeData(wb, "Overview", data.frame(c(title, description)),
    startRow = 1, colNames = FALSE
  )
  addStyle(wb, "Overview", createStyle(textDecoration = "bold"), rows = 1, cols = 1)
  writeData(wb, "Overview", index, startRow = 4)
  addStyle(wb, "Overview", createStyle(textDecoration = "bold", fgFill = "#DCE6F1"),
    rows = 4, cols = seq_len(ncol(index)), gridExpand = TRUE
  )
  setColWidths(wb, "Overview", cols = 1:3, widths = "auto")
  setColWidths(wb, "Overview", cols = 4, widths = 90)
  worksheetOrder(wb) <- c(length(present) + 1L, seq_along(present))
  invisible(index)
}

# Every stage workbook writes its sheets the same way: bold header on the pale
# fill, frozen top row, columns sized to their contents. Four scripts carried
# byte-identical copies of this and a fifth differed by one line. The remove
# guard is that line, and it lets a script that loads an existing workbook
# re-run without colliding with the sheets it wrote last time.
write_sheet <- function(wb, name, data) {
  if (name %in% names(wb)) removeWorksheet(wb, name)
  addWorksheet(wb, name)
  writeData(wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
  )
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}
