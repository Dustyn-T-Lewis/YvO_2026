# Print-scale override for the oversized main composite.
#
# A composite authored wider than the print column is shrunk linearly when it
# is placed, so its text has to be enlarged by the same factor to survive. The
# factor was hard-coded at 380/178 while F04, the only caller left since F05
# retired, is authored at 460 mm: panel text was landing at 83% of its intended
# size. It now follows the canvas.
#
# WARNING: this script MUTATES style.R's exported globals
# (FIG_TITLE_SIZE, FIG_SUBTITLE_SIZE, FIG_STRIP_SIZE, FIG_AXIS_TEXT,
#  FIG_LEGEND_TITLE, FIG_LEGEND_TEXT, FIG_THEME, BASE_*) and never restores
# them. Anything sourced afterwards at a different width has to re-source
# style.R first; 04_Figures/F04/a_script/F04.R does exactly that.
#
# Source order: style.R first, this second, then the panel scripts.

PRINT_CANVAS_MM <- if (exists("PRINT_CANVAS_MM")) PRINT_CANVAS_MM else 460
PRINT_SCALE <- PRINT_CANVAS_MM / 178

# The BASE_* annotation sizes keep the older, smaller factor. They are not free
# text: panel A's label boxes are fitted to them, and enlarging them by the full
# canvas factor pushes the pathway bars over the central scatter. The FIG_*
# theme sizes carry no layout constraint and do follow the canvas, which is
# where the unreadable axis text was.
LABEL_SCALE <- 380 / 178

FIG_TITLE_SIZE    <- round(FIG_TITLE_SIZE * PRINT_SCALE, 1)
FIG_SUBTITLE_SIZE <- round(FIG_SUBTITLE_SIZE * PRINT_SCALE, 1)
FIG_STRIP_SIZE    <- round(FIG_STRIP_SIZE * PRINT_SCALE, 1)
FIG_AXIS_TEXT     <- round(FIG_AXIS_TEXT * PRINT_SCALE, 1)
FIG_LEGEND_TITLE  <- round(FIG_LEGEND_TITLE * PRINT_SCALE, 1)
FIG_LEGEND_TEXT   <- round(FIG_LEGEND_TEXT * PRINT_SCALE, 1)

# drop 2pt to avoid overfilling
BASE_PATHWAY  <- round(BASE_PATHWAY * LABEL_SCALE, 1) - 2
BASE_GENE     <- round(BASE_GENE * LABEL_SCALE, 1) - 2
BASE_STAT     <- round(BASE_STAT * LABEL_SCALE, 1) - 2
BASE_QUADRANT <- round(BASE_QUADRANT * LABEL_SCALE, 1) - 2
BASE_COUNT    <- round(BASE_COUNT * LABEL_SCALE, 1) - 2

FIG_THEME <- theme_bw(base_size = round(6 * PRINT_SCALE, 1), base_family = "Helvetica") +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      color = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag           = element_text(face = "bold", size = FIG_TITLE_SIZE),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(t = 0)),
    axis.title.y       = element_text(face = "bold", size = FIG_AXIS_TEXT,
                                      margin = margin(r = -1)),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(face = "bold", size = FIG_LEGEND_TEXT,
                                      color = "grey15"),
    legend.key.size    = unit(2.5 * PRINT_SCALE, "mm"),
    panel.grid.minor   = element_blank()
  )
