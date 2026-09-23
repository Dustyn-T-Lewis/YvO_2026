#!/usr/bin/env Rscript
# The WGCNA module palette and the helpers that keep type readable on it.
#
# Its own file, with no side effects, because abstract_panels/_common.R needs
# it but cannot source style.R: that executes devices.R and defines size
# globals print_scale_apply.R mutates in place. style.R sources this, so every
# F05/F06 script gets it for free.

# The nine WGCNA modules, drawn. The module is named by its colour, so the
# colour is an identifier and has to be the same in every figure that shows a
# module -- which is why this lives here and not in a panel script. These are
# the raw WGCNA names moved far enough to survive a printer: pure yellow and
# pure green have almost no ink, and pure black leaves no headroom for a label.
# The hue is preserved, so a module still reads as the colour it is named for.
MODULE_FILL <- c(
  turquoise = "#3FBFB0", blue = "#3779AC", brown = "#8C6239",
  yellow = "#E3B505", green = "#4CA64C", black = "#3A3A3A",
  red = "#CE4948", pink = "#E48FA8", magenta = "#B5478F",
  grey = "#BFBFBF"
)

# Short tags for places a bio label will not fit.
MODULE_TAG <- c(
  turquoise = "Fatty acid ox.", blue = "Chaperones", brown = "Proteasome",
  yellow = "Sarcomere", green = "OxPhos", black = "Glycolysis",
  red = "Ribosome", pink = "ATP synthase", magenta = "40S ribosome"
)

# Returning NA for an unknown module would draw a transparent bar and say
# nothing, so an unknown module is an error.
module_fill <- function(module) {
  key <- sub("^ME", "", module)
  hex <- unname(MODULE_FILL[key])
  bad <- is.na(hex) & !is.na(key)
  if (any(bad)) {
    stop("no MODULE_FILL entry for: ", paste(unique(key[bad]), collapse = ", "))
  }
  hex
}

# WCAG AA for body text, and the dark ink on_fill weighs against white.
WCAG_AA <- 4.5
DARK_INK <- "grey10"

relative_luminance <- function(hex) {
  v <- col2rgb(hex) / 255
  f <- ifelse(v <= 0.03928, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  as.numeric(0.2126 * f[1, ] + 0.7152 * f[2, ] + 0.0722 * f[3, ])
}

# Ink that stays readable on a given fill. Compare both candidates and take
# whichever wins on contrast, rather than testing luminance against a cut. A
# cut is what put white type on green at 3.06:1 when dark type gives 5.69:1 --
# green sits within a hundredth of any plausible threshold, which is the same
# trap the hand-kept list of pale modules fell into.
on_fill <- function(hex) {
  lum <- relative_luminance(hex)
  cr_dark <- (lum + 0.05) / (relative_luminance(DARK_INK) + 0.05)
  cr_light <- 1.05 / (lum + 0.05)
  ifelse(cr_dark >= cr_light, DARK_INK, "white")
}

# A fill darkened until type drawn in it clears a contrast ratio on white,
# keeping the hue and compensating saturation so the module stays recognisable.
# Raw yellow as text is invisible; this is what makes a module colour usable as
# a label rather than only as a block.
text_safe <- function(hex) {
  vapply(hex, function(h) {
    hsv_in <- rgb2hsv(col2rgb(h))
    for (k in seq(1, 0, by = -0.005)) {
      sat <- min(1, hsv_in[2] * (2 - k))
      out <- hsv(hsv_in[1], sat, hsv_in[3] * k)
      if (1.05 / (relative_luminance(out) + 0.05) >= WCAG_AA) {
        return(out)
      }
    }
    "#000000"
  }, character(1), USE.NAMES = FALSE)
}
