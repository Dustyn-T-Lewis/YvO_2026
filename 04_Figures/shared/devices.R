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

# Writes <key>_Figure.pdf, the supplementary figure over its legend from
# supp_legends.txt: 178 mm column, 16 mm margins, 8 pt legend on A4. A figure
# too large for that column gets a proportionally larger A4-shaped page instead
# of being redrawn smaller, so the figure itself is drawn exactly as saved.
# `fig` is one plot or grob, or the path of a PDF that already holds the pages.
# A multi-page PDF, or a figure that shrinking would leave under half the
# column wide, gets the legend on a page of its own, then the figure at its size.
caption_supp <- function(fig, key, width_mm, height_mm, dir) {
  legends <- strsplit(paste(readLines("04_Figures/shared/supp_legends.txt"), collapse = "\n"), "\n\n")[[1]]
  legend <- legends[startsWith(legends, paste0(key, " Figure."))]
  stopifnot("no legend for this figure in supp_legends.txt" = length(legend) == 1)

  legend_grob <- function(k) {
    gridtext::textbox_grob(legend,
      x = grid::unit(16 * k, "mm"), y = grid::unit(1, "npc"), hjust = 0, vjust = 1,
      width = grid::unit(178 * k, "mm"),
      gp = grid::gpar(fontsize = 8 * k, fontfamily = "Helvetica", lineheight = 1.25)
    )
  }
  pdf(NULL)
  legend_h <- grid::convertHeight(grid::grobHeight(legend_grob(1)), "mm", valueOnly = TRUE)
  dev.off()

  draw <- function(p) if (inherits(p, "grob")) grid::grid.draw(p) else print(p, newpage = FALSE)
  k <- max(1, width_mm / 178, height_mm / (297 - 32 - legend_h - 4))
  dev_fun <- get_pdf_device()

  if (is.character(fig) || (k > 1 && width_mm / k < 89)) {
    legend_page <- tempfile(fileext = ".pdf")
    dev_fun(legend_page, width = 210 / 25.4, height = 297 / 25.4)
    grid::pushViewport(grid::viewport(y = grid::unit(1, "npc") - grid::unit(16, "mm"), just = "top"))
    grid::grid.draw(legend_grob(1))
    dev.off()
    if (!is.character(fig)) {
      path <- tempfile(fileext = ".pdf")
      dev_fun(path, width = width_mm / 25.4, height = height_mm / 25.4)
      draw(fig)
      dev.off()
      fig <- path
    }
    qpdf::pdf_combine(c(legend_page, fig), file.path(dir, paste0(key, "_Figure.pdf")))
    return(invisible())
  }

  dev_fun(file.path(dir, paste0(key, "_Figure.pdf")), width = 210 * k / 25.4, height = 297 * k / 25.4)
  top <- grid::unit(1, "npc") - grid::unit(16 * k, "mm")
  grid::pushViewport(grid::viewport(y = top, width = grid::unit(width_mm, "mm"),
                                    height = grid::unit(height_mm, "mm"), just = "top"))
  draw(fig)
  grid::popViewport()
  grid::pushViewport(grid::viewport(y = top - grid::unit(height_mm + 4 * k, "mm"), just = "top"))
  grid::grid.draw(legend_grob(k))
  dev.off()
  invisible()
}
