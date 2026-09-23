# Device selection, shared by the stage report scripts and the figure tree.
#
# Kept separate from style.R so a stage report can pick a working PDF device
# without also loading the figure palettes and theme.

# A cairo build without its X11 module signals "failed to load cairo DLL" as a
# warning, not an error, and then writes nothing. Catching only errors let that
# broken device through and every PDF save failed silently, which is how the
# F03 significance heatmaps shipped as PNG alone. A warning counts as failure.
probe_device <- function(expr) {
  tryCatch(
    {
      expr()
      dev.off()
      TRUE
    },
    error = function(e) FALSE,
    warning = function(w) {
      while (length(dev.list())) dev.off()
      FALSE
    }
  )
}

get_pdf_device <- function() {
  # cairo_pdf > quartz > base pdf.
  #
  # quartz sits ahead of the base device because its font repertoire is far
  # wider. Measured 2026-09-16 on the same label: quartz keeps the Greek and
  # the dashes the figures are full of, and loses three characters; the base
  # device loses eleven, including every Pi, Delta and rho. Preferring the base
  # device would turn six holes into about sixty-five.
  if (probe_device(function() cairo_pdf(tempfile()))) {
    return(cairo_pdf)
  }

  if (probe_device(function() quartz(type = "pdf", file = tempfile(fileext = ".pdf")))) {
    return(function(filename, width, height, ...) {
      quartz(file = filename, type = "pdf", width = width, height = height)
    })
  }

  "pdf"
}

get_raster_pdf_device <- function() {
  # For composites assembled from pre-rendered panel PNGs, which the quartz
  # device draws upside down: Figure 5, S7 Figure and S8 Figure all shipped
  # mirrored while their PNG twins were correct, and nothing warned. The base
  # device places rasters the right way up.
  #
  # Only the three raster composites use it, because of the repertoire gap
  # described above. They carry no text beyond panel letters and two titles,
  # so losing the wider glyph set costs them nothing.
  if (probe_device(function() cairo_pdf(tempfile()))) {
    return(cairo_pdf)
  }

  grDevices::pdf
}

# The base device names its fonts rather than carrying them, and a submission
# system reads a named-only font as a missing one. Ghostscript rewrites the
# finished file with the fonts attached, which leaves how the figure was drawn
# alone -- the point, since the raster composites are on that device precisely
# to keep their panels the right way up.
#
# grDevices::embedFonts() is not enough on its own: Helvetica is one of the
# PDF base-14, which Ghostscript keeps on its NeverEmbed list and so leaves
# referenced by name. Clearing that list is the whole fix. Doing it this way
# rather than through -dPDFSETTINGS=/prepress, which also embeds, is
# deliberate: prepress downsampled panel A from 512 to 300 ppi.
embed_pdf_fonts <- function(path) {
  gs <- tools::find_gs_cmd()
  if (!nzchar(gs)) {
    warning(
      "Ghostscript not found; ", basename(path), " keeps unembedded fonts",
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  staged <- paste0(path, ".embed")
  status <- system2(gs, c(
    "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
    "-dEmbedAllFonts=true", "-dSubsetFonts=true",
    paste0("-sOutputFile=", shQuote(staged)),
    "-c", shQuote("<</NeverEmbed []>> setdistillerparams"),
    "-f", shQuote(path)
  ), stdout = FALSE, stderr = FALSE)
  if (status != 0 || !file.exists(staged) || file.size(staged) == 0) {
    unlink(staged)
    warning(
      "Ghostscript could not embed fonts in ", basename(path),
      call. = FALSE
    )
    return(invisible(FALSE))
  }
  file.rename(staged, path)
  invisible(TRUE)
}

# Under Rscript, an implicit device open writes Rplots.pdf into the working
# directory. Redirect that device to a null sink; every real output goes
# through an explicit device, so nothing is lost.
suppress_rplots_device <- function() {
  options(device = function(...) grDevices::pdf(file = nullfile(), ...))
}

suppress_rplots_device()

# Multi-page reports open a device, print into it and close it, rather than
# going through ggsave(device = ). Resolve the device to a function so the
# base-pdf fallback works the same way as the other two.
open_pdf <- function(path, width, height) {
  dev <- get_pdf_device()
  if (is.character(dev)) dev <- match.fun(dev)
  dev(path, width = width, height = height)
}
