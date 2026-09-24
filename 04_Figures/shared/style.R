# style.R — palettes, themes, sizing, helpers for all figure scripts
#
# Everything that is a repo-wide constant lives in tree_config.R.

source("04_Figures/shared/tree_config.R")
source("04_Figures/shared/module_palette.R")

pacman::p_load(ggplot2, scales, grid)

source("04_Figures/shared/devices.R")

AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")

DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  Young_Pre  = scales::alpha("#4393C3", 0.5),
  Young_Post = "#4393C3",
  Old_Pre    = scales::alpha("#D6604D", 0.5),
  Old_Post   = "#D6604D"
)

SHAPE_TP <- c(Pre = 16, Post = 17)

CONTRAST_COLORS <- c(
  Aging          = "#4CAF50",
  Training_Young = "#E05A4E",
  Training_Old   = "#5DA5DA",
  Interaction    = "#9B7FBF"
)

DB_COLORS <- c(
  Hallmark = "#AA336A", KEGG = "#E65100", Reactome = "#1565C0",
  "GO:BP" = "#00796B",
  "GO Slim" = "#8D6E63", BioCarta = "#795548", PID = "#455A64"
)

PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#377EB8", MNAR = "#E41A1C")

PCA_COLORS <- c(
  Young_Pre  = "#93C4DE",
  Young_Post = "#2166AC",
  Old_Pre    = "#F4A582",
  Old_Post   = "#B2182B"
)

PCA_SHAPES <- c(Young_Pre = 16, Young_Post = 17, Old_Pre = 16, Old_Post = 17)

SIG_COLORS_F2 <- c(
  "Interaction"    = "#7B5EA7",
  "Sig Both"       = "#2E7D32",
  "Sig Young only" = "#E05A4E",
  "Sig Old only"   = "#5DA5DA",
  "NS"             = "grey70"
)

SIG_LABEL_FILL_F2 <- c(
  "Interaction"    = scales::alpha("#7B5EA7", 0.75),
  "Sig Both"       = scales::alpha("#2E7D32", 0.75),
  "Sig Young only" = scales::alpha("#E05A4E", 0.75),
  "Sig Old only"   = scales::alpha("#5DA5DA", 0.75),
  "NS"             = scales::alpha("grey70", 0.75)
)
SIG_LABEL_TEXT_F2 <- setNames(rep("white", 5), names(SIG_LABEL_FILL_F2))

SIG_COLORS_F3 <- c(
  "Sig Both"           = "#2E7D32",
  "Sig Aging only"     = "#E05A4E",
  "Sig Training only"  = "#5DA5DA",
  "NS"                 = "grey70"
)

SIG_LABEL_FILL_F3 <- c(
  "Sig Both"           = scales::alpha("#2E7D32", 0.75),
  "Sig Aging only"     = scales::alpha("#E05A4E", 0.75),
  "Sig Training only"  = scales::alpha("#5DA5DA", 0.75),
  "NS"                 = scales::alpha("grey70", 0.75)
)
SIG_LABEL_TEXT_F3 <- setNames(rep("white", 4), names(SIG_LABEL_FILL_F3))

ORA_QUAD_COLORS_F2 <- c(
  "Concordant Up"               = "#E57373",
  "Concordant Down"             = "#64B5F6",
  "Discordant (Y Up / O Down)"  = "#FFB74D",
  "Discordant (Y Down / O Up)"  = "#81C784"
)

ORA_QUAD_COLORS_F3 <- c(
  "Reversed (Aging Up / Training Down)"  = "#E57373",
  "Reversed (Aging Down / Training Up)"  = "#64B5F6",
  "Exacerbated Up"                       = "#FFB74D",
  "Exacerbated Down"                     = "#81C784"
)

CLUSTER_COLORS <- c(
  C1 = "#E74C3C", C2 = "#3498DB", C3 = "#2ECC71",
  C4 = "#F39C12", C5 = "#9B59B6", C6 = "#1ABC9C",
  C7 = "#E67E22", C8 = "#34495E", C9 = "#D35400",
  C10 = "#7F8C8D"
)

THEME_COLORS <- c(
  "Mitochondrial & Energy Metabolism" = "#E57373",
  "Muscle Structure & Myogenesis" = "#64B5F6",
  "Proteostasis & Stress Response" = "#81C784",
  "Cytoskeletal & Cell Division" = "#CE93D8",
  "Immune & Complement" = "#FFB74D",
  "ECM & Tissue Remodeling" = "#F05292",
  "Metabolic & Redox Regulation" = "#FFD54F",
  "Intracellular Transport & Signaling" = "#4DB6AC"
)

PANEL_MD <- 178 # J Physiol double-column width

# base annotation sizes (mm) at PANEL_MD; scale_text() adjusts for other widths
BASE_PATHWAY <- 2.8 # ~8pt pathway labels
BASE_GENE <- 2.5 # ~7pt gene labels
BASE_STAT <- 2.5 # ~7pt stat annotations
BASE_QUADRANT <- 2.8 # ~8pt quadrant labels
BASE_COUNT <- 2.5 # ~7pt bar counts
BASE_TAG <- 8 # panel tags

# Linear, not sqrt. A composite drawn on a wide canvas is shrunk linearly when
# it is placed at ref_width, so a square-root enlargement undoes only part of
# it: on F04's 460 mm canvas a 7 pt annotation used to print at 4.4 pt.
scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * (panel_width_mm / ref_width)
}

# strip decorations before embedding in composite (call after standalone save)
strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}
# text hierarchy — J Physiol spec
FIG_TITLE_SIZE <- 7
FIG_SUBTITLE_SIZE <- 6
FIG_STRIP_SIZE <- 5
FIG_AXIS_TEXT <- 5
FIG_LEGEND_TITLE <- 5
FIG_LEGEND_TEXT <- 6

# Tag, title and subtitle sizes for a composite.
#
# Keyed off canvas WIDTH against the width the figure will be printed at, not
# canvas height. Height was the wrong input: it made a short wide composite and
# a tall narrow one disagree on title size even though both print at 178 mm.
#
# Most figures are authored at their print width, so `print_w_mm` defaults to
# the canvas and the sizes come out at the constants below. F04 and F05 are
# authored oversized and shrink at print, so they pass the print width
# explicitly and their type scales up to compensate.
composite_text_sizes <- function(canvas_w_mm, print_w_mm = canvas_w_mm) {
  k <- canvas_w_mm / print_w_mm
  list(
    title    = round(FIG_TITLE_SIZE * k, 1),
    subtitle = round(FIG_SUBTITLE_SIZE * k, 1),
    tag      = round(BASE_TAG * k, 1)
  )
}

FIG_THEME <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    plot.title = element_text(
      face = "bold", size = FIG_TITLE_SIZE,
      margin = margin(b = 1)
    ),
    plot.subtitle = element_text(
      face = "bold.italic", size = FIG_SUBTITLE_SIZE,
      color = "grey30", margin = margin(t = 0, b = 2)
    ),
    plot.tag = element_text(face = "bold", size = BASE_TAG),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x = element_text(
      face = "bold", size = 5,
      margin = margin(t = 0)
    ),
    axis.title.y = element_text(
      face = "bold", size = 5,
      margin = margin(r = -1)
    ),
    axis.text = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title = element_text(
      face = "bold", size = FIG_LEGEND_TITLE,
      color = "grey20"
    ),
    legend.text = element_text(
      face = "bold", size = FIG_LEGEND_TEXT,
      color = "grey15"
    ),
    legend.key.size = unit(2.5, "mm"),
    panel.grid.minor = element_blank()
  )

fmt_p <- function(p) {
  if (p < 0.001) {
    return("p < 0.001")
  }
  if (p < 0.01) {
    return(sprintf("p = %.3f", p))
  }
  sprintf("p = %.2f", p)
}

# bolditalic(), not the bold(italic()) nesting these used to emit: the quartz
# PDF device and the PNG device both drop the outer bold from the nested form,
# so every p the two helpers meant to emphasise printed at the same weight as
# the ones they did not.
fmt_p_plot <- function(p, threshold = 0.05) {
  label <- fmt_p(p)
  if (p < threshold) sprintf('bolditalic("%s")', label) else sprintf('italic("%s")', label)
}

fmt_anova_sub <- function(age_p, time_p, int_p, threshold = 0.05) {
  wrap <- function(label, p) {
    txt <- sprintf("%s %s", label, fmt_p(p))
    if (p < threshold) sprintf('bolditalic("%s")', txt) else sprintf('italic("%s")', txt)
  }
  paste('italic("RM-ANOVA:")', wrap("Age", age_p), wrap("Time", time_p),
    wrap("Int.", int_p),
    sep = " ~~ "
  )
}

# Bonett & Wright 2000 — Fisher z CI for r (k = number of covariates)
fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) {
    return(c(lo = NA_real_, hi = NA_real_))
  }
  z <- atanh(r)
  se <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

# Returns a one-row data frame so summarise() unpacks it into ci_lo/ci_hi.
# Indexing a named vector twice, as this was called before, ran the bootstrap
# twice and took each bound from a different draw.
#
# Draws from the caller's RNG rather than seeding itself: both callers seed at
# their top, and seeding here would shift the published intervals. A new caller
# must seed before calling, or its intervals move every run.
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs <- unname(quantile(meds, c((1 - conf) / 2, (1 + conf) / 2)))
  data.frame(ci_lo = qs[1], ci_hi = qs[2])
}

.DB_PREFIXES <- c(
  "^HALLMARK_", "^GOSLIM_", "^GOBP_", "^GOCC_", "^GOMF_",
  "^REACTOME_", "^KEGG_MEDICUS_", "^KEGG_"
)

.SCI_CAPS <- c(
  "Mtorc1" = "mTORC1", "Myc " = "MYC ", "E2f " = "E2F ",
  "Dna " = "DNA ", "Rna " = "RNA ", "Tnfa " = "TNFa ",
  "Uv " = "UV ", "G2m " = "G2M ", "Il6 " = "IL6 ",
  "Il2 " = "IL2 ", "Kras " = "KRAS ", "P53 " = "p53 ",
  "Tgf " = "TGF ", "Nfkb" = "NF-kB", "Atp " = "ATP ",
  "Nadh " = "NADH ",
  "Oxidative Phosphorylation" = "OXPHOS",
  "External Encapsulating Structure Or.*" = "Extracellular Matrix Organization",
  "Enzyme Linked Receptor Protein Signaling.*" = "Receptor Protein Signaling"
)

clean_pathway_name <- function(name) {
  out <- name
  for (pfx in .DB_PREFIXES) out <- stringr::str_remove(out, pfx)
  out <- stringr::str_replace_all(out, "_", " ")
  out <- stringr::str_to_title(out)
  for (i in seq_along(.SCI_CAPS)) {
    out <- stringr::str_replace(out, names(.SCI_CAPS)[i], .SCI_CAPS[i])
  }
  out
}

CTR_SHORT <- c(
  Aging          = "Aging",
  Training_Young = "Tr.(Y)",
  Training_Old   = "Tr.(O)",
  Interaction    = "Inter."
)
CTR_FACET <- CTR_SHORT
CONTRAST_ORDER <- c("Aging", "Training_Young", "Training_Old", "Interaction")
