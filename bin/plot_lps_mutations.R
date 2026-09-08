#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# plot_lps_mutations.R
#
# Publication figures for the Pasteurella multocida LPS typing pipelines
# (Illumina and ONT).  For each LPS type present in a run it draws the same
# information as the lollipop plots in the HTML report -- one lollipop per
# distinct mutation, height = number of genomes carrying it -- above a
# directional gene-arrow map of the reference locus, and writes it as PDF/SVG.
#
# Run by hand against a finished results directory.  It is not part of the
# Nextflow workflow and needs no container.
#
# Dependencies (CRAN):
#   install.packages(c("ggplot2", "dplyr", "scales", "svglite"))
# or re-run this script with --install-deps to install whatever is missing.
#
# Example:
#   Rscript bin/plot_lps_mutations.R \
#     --report-dir results/10_report \
#     --lps-db-dir databases/LPS \
#     --outdir figures
# ---------------------------------------------------------------------------

REQUIRED_PKGS <- c("ggplot2", "dplyr", "scales", "svglite")

# --------------------------------------------------------------------------- #
# Command line
# --------------------------------------------------------------------------- #

# options taking no value
FLAG_OPTS <- c("help", "install-deps", "write-table")

VALUE_OPTS <- c("report-dir", "subtype-report", "lps-db-dir", "outdir", "types",
                "samples", "samples-file", "exclude-samples", "palette",
                "gene-colors", "min-genomes", "legend", "width", "height",
                "font", "prefix", "formats")

USAGE <- "
Usage: plot_lps_mutations.R --lps-db-dir DIR [options]

Input
  --report-dir DIR        directory holding 10_*_subtype_report.tsv (default '.').
                          Comma-separated list accepted.
  --subtype-report FILE   explicit report path(s), comma-separated. Overrides
                          --report-dir.
  --lps-db-dir DIR        LPS database: reference_LPS.txt, *.gb, gene_colors.tsv.
                          Required.

Selection
  --types L3,L6           only these LPS types (default: all)
  --samples PM1,PM2       only these genomes (default: all)
  --samples-file FILE     one sample id per line
  --exclude-samples PM3   drop these genomes
  --min-genomes N         drop mutations seen in fewer than N genomes (default 1)

Appearance
  --palette db|okabe-ito  gene colours (default db = gene_colors.tsv)
  --gene-colors FILE      alternative GENE<tab>HEX table
  --legend all|unlabelled|none
                          gene legend with coordinates (default none; every gene
                          is already labelled on the track). 'unlabelled' lists
                          only genes too narrow to label inside their arrow
  --font NAME             font family (default Helvetica)
  --width N --height N    figure size in inches (default: scaled to the locus)

Output
  --outdir DIR            output directory (default lps_mutation_figures)
  --prefix STR            filename prefix (default lps_mutations_)
  --formats pdf,svg       output formats, any of pdf,svg,png (default pdf,svg)
  --write-table           also write the plotted variants as a TSV

Other
  --install-deps          install missing R packages from CRAN, then continue
  --help                  show this message
"

parse_args <- function(argv) {
  known <- c(FLAG_OPTS, VALUE_OPTS)
  opts <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("unexpected argument: ", a, call. = FALSE)
    key <- sub("^--", "", a)
    val <- NULL
    if (grepl("=", key, fixed = TRUE)) {  # --key=value form
      val <- sub("^[^=]*=", "", key)
      key <- sub("=.*$", "", key)
    }
    if (!key %in% known) stop("unknown option --", key, "\n", USAGE, call. = FALSE)
    if (key %in% FLAG_OPTS) {
      opts[[key]] <- TRUE
    } else {
      if (is.null(val)) {
        if (i == length(argv)) stop("missing value for --", key, call. = FALSE)
        i <- i + 1L
        val <- argv[[i]]
      }
      opts[[key]] <- val
    }
    i <- i + 1L
  }
  opts
}

args_list <- parse_args(commandArgs(trailingOnly = TRUE))

opt <- function(name, default = NULL) {
  v <- args_list[[name]]
  if (is.null(v)) default else v
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

if (isTRUE(opt("help"))) { cat(USAGE); quit(status = 0) }

# --------------------------------------------------------------------------- #
# Dependencies
# --------------------------------------------------------------------------- #

ensure_packages <- function(pkgs, install = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible(NULL))
  install_line <- sprintf('install.packages(c(%s))',
                          paste(sprintf('"%s"', missing), collapse = ", "))
  if (!install) {
    stop("missing R package(s): ", paste(missing, collapse = ", "),
         "\nInstall them with:\n  ", install_line,
         "\nor re-run with --install-deps", call. = FALSE)
  }
  message("installing: ", paste(missing, collapse = ", "))
  utils::install.packages(missing, repos = "https://cloud.r-project.org")
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still)) {
    stop("could not install: ", paste(still, collapse = ", "),
         "\nInstall manually with:\n  ", install_line, call. = FALSE)
  }
  invisible(NULL)
}

ensure_packages(REQUIRED_PKGS, install = isTRUE(opt("install-deps")))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

NA_VALUES <- c("", "NA", "na", "N/A", ".", "-")

is_blank <- function(x) is.na(x) | trimws(x) %in% NA_VALUES

# Fallback gene-track palette, kept identical to generate_lps_report.py so the
# report and these figures colour unknown genes the same way.
FALLBACK_PALETTE <- c(
  "#D85A30", "#378ADD", "#7F77DD", "#639922", "#BA7517",
  "#E24B4A", "#1D9E75", "#D4537E", "#888780", "#4682B4",
  "#8B4513", "#DA70D6"
)

OKABE_ITO <- c("#0072B2", "#009E73", "#E69F00", "#D55E00", "#CC79A7",
               "#56B4E9", "#F0E442", "#999999", "#7F3C8D", "#3969AC",
               "#11A579", "#E68310")

INTERGENIC_COLOUR <- "#5F5E5A"

read_table_tsv <- function(path, header = TRUE) {
  utils::read.delim(path, sep = "\t", quote = "", comment.char = "",
                    colClasses = "character", check.names = FALSE,
                    header = header, na.strings = character(0))
}

# Shorten a VARTYPE for plotting: 'snp S92*' -> 'S92*', '19bp deletion' -> '19bp del'.
# The database is inconsistent about the space after 'snp' (e.g. 'snpY226*').
short_vartype <- function(v) {
  v <- trimws(ifelse(is.na(v), "", v))
  v <- sub("^snp\\s*", "", v, perl = TRUE, ignore.case = TRUE)
  v <- gsub("\\binsertion\\b", "ins", v, perl = TRUE, ignore.case = TRUE)
  v <- gsub("\\bdeletion\\b", "del", v, perl = TRUE, ignore.case = TRUE)
  trimws(v)
}

# Black or white text, whichever has more contrast against `hex`.
contrast_text <- function(hex) {
  m <- grDevices::col2rgb(hex) / 255
  lin <- ifelse(m <= 0.03928, m / 12.92, ((m + 0.055) / 1.055)^2.4)
  lum <- 0.2126 * lin[1, ] + 0.7152 * lin[2, ] + 0.0722 * lin[3, ]
  ifelse(lum > 0.45, "grey10", "white")
}

# Rendered width of each string, in inches, measured on a null device.
measure_widths <- function(txt, pointsize, font = 1L, family = "") {
  if (!length(txt)) return(numeric(0))
  grDevices::pdf(NULL, width = 12, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  ok <- tryCatch({ graphics::par(family = family); TRUE },
                 error = function(e) FALSE, warning = function(e) FALSE)
  if (!ok) graphics::par(family = "")
  graphics::par(ps = pointsize, font = font)
  graphics::strwidth(txt, units = "inches")
}

# --------------------------------------------------------------------------- #
# LPS database
# --------------------------------------------------------------------------- #

LOCUS_LEN_RE <- "^LOCUS\\s+\\S+\\s+(\\d+)\\s+bp"
GENE_LOC_RE  <- "^\\s{1,8}gene\\s+(complement\\()?\\s*<?(\\d+)\\.\\.>?(\\d+)\\)?\\s*$"
GENE_NAME_RE <- "/gene=\"([^\"]+)\""

# Return list(length, genes) from a GenBank file.  Only top-level `gene`
# features are used; start/end are 1-based inclusive, strand is '+' or '-'
# (taken from the complement() wrapper).  Mirrors parse_genbank_genes() in
# generate_lps_report.py so both renderers agree on gene direction.
parse_genbank_genes <- function(path) {
  lines <- readLines(path, warn = FALSE)
  locus_len <- 0L
  rows <- list()
  pending <- NULL
  for (ln in lines) {
    if (locus_len == 0L) {
      m <- regmatches(ln, regexec(LOCUS_LEN_RE, ln, perl = TRUE, ignore.case = TRUE))[[1]]
      if (length(m)) locus_len <- as.integer(m[2])
    }
    m <- regmatches(ln, regexec(GENE_LOC_RE, ln, perl = TRUE))[[1]]
    if (length(m)) {
      pending <- list(start = as.integer(m[3]), end = as.integer(m[4]),
                      strand = if (nzchar(m[2])) "-" else "+")
      next
    }
    if (!is.null(pending)) {
      nm <- regmatches(ln, regexec(GENE_NAME_RE, ln, perl = TRUE))[[1]]
      if (length(nm)) {
        rows[[length(rows) + 1L]] <- data.frame(
          name = trimws(nm[2]), start = pending$start, end = pending$end,
          strand = pending$strand, stringsAsFactors = FALSE)
        pending <- NULL
      }
    }
  }
  if (!length(rows)) return(list(length = locus_len, genes = NULL))
  # `gene` and `CDS` both carry /gene=, so the same gene can appear twice
  genes <- do.call(rbind, rows)
  genes <- genes[!duplicated(genes[c("name", "start", "end")]), , drop = FALSE]
  genes <- genes[order(genes$start), , drop = FALSE]
  rownames(genes) <- NULL
  if (!all(genes$strand %in% c("+", "-"))) {
    stop("could not determine strand for every gene in ", path, call. = FALSE)
  }
  if (locus_len == 0L) locus_len <- max(genes$end)
  list(length = locus_len, genes = genes)
}

read_reference_map <- function(db_dir) {
  path <- file.path(db_dir, "reference_LPS.txt")
  if (!file.exists(path)) stop("no reference_LPS.txt in ", db_dir, call. = FALSE)
  tbl <- read_table_tsv(path, header = FALSE)  # TYPE / .gb / .fasta, no header
  if (ncol(tbl) < 2) stop("reference_LPS.txt needs at least 2 columns", call. = FALSE)
  keep <- nzchar(trimws(tbl[[1]]))
  stats::setNames(trimws(tbl[[2]])[keep], trimws(tbl[[1]])[keep])
}

read_gene_colors <- function(path) {
  if (is.null(path) || !file.exists(path)) return(character(0))
  lines <- readLines(path, warn = FALSE)
  out <- character(0)
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    if (length(parts) < 2 || !nzchar(trimws(parts[1]))) next
    if (i == 1L && tolower(trimws(parts[2])) %in% c("hex", "color", "colour")) next
    out[trimws(parts[1])] <- trimws(parts[2])
  }
  out
}

# Colour every gene in locus order: database colour where available, otherwise
# the deterministic fallback palette.
color_for_genes <- function(gene_names, gene_colors, palette) {
  if (identical(palette, "okabe-ito")) {
    return(stats::setNames(OKABE_ITO[(seq_along(gene_names) - 1L) %% length(OKABE_ITO) + 1L],
                           gene_names))
  }
  out <- character(0)
  pi <- 0L
  for (nm in gene_names) {
    if (nm %in% names(gene_colors)) {
      out[nm] <- gene_colors[[nm]]
    } else {
      out[nm] <- FALLBACK_PALETTE[pi %% length(FALLBACK_PALETTE) + 1L]
      pi <- pi + 1L
    }
  }
  out
}

# --------------------------------------------------------------------------- #
# Subtype report
# --------------------------------------------------------------------------- #

NEEDED_COLS <- c("SAMPLE", "TYPE", "VARTYPE", "CHROM", "POS", "REF", "ALT", "GENE")

read_subtype_reports <- function(paths) {
  parts <- lapply(paths, function(p) {
    df <- read_table_tsv(p)
    missing <- setdiff(NEEDED_COLS, names(df))
    if (length(missing)) {
      stop(p, " is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    df[NEEDED_COLS]
  })
  df <- do.call(rbind, parts)
  # a genome sequenced on two platforms must not be counted twice
  df[!duplicated(df[c("SAMPLE", "CHROM", "POS", "REF", "ALT")]), , drop = FALSE]
}

discover_reports <- function(dirs) {
  found <- unlist(lapply(dirs, function(d) {
    Sys.glob(file.path(d, "10_*_subtype_report.tsv"))
  }))
  if (!length(found)) {
    stop("no 10_*_subtype_report.tsv found in: ", paste(dirs, collapse = ", "),
         call. = FALSE)
  }
  found
}

# --------------------------------------------------------------------------- #
# Geometry
# --------------------------------------------------------------------------- #

# Five-point gene arrow, drawn counter-clockwise and flipped for the minus strand.
arrow_poly <- function(x1, x2, y, h, strand, head) {
  hl <- min(head, 0.5 * (x2 - x1))
  if (identical(strand, "-")) {
    data.frame(x = c(x2, x1 + hl, x1, x1 + hl, x2),
               y = c(y - h, y - h, y, y + h, y + h))
  } else {
    data.frame(x = c(x1, x2 - hl, x2, x2 - hl, x1),
               y = c(y - h, y - h, y, y + h, y + h))
  }
}

# Place each mutation label just above its own marker, lifting it by `dy` only
# where it would collide with a label already placed or ride over a taller
# neighbouring lollipop.  Left to right and deterministic, unlike ggrepel: an
# isolated label keeps a short leader line, and only clusters stack up.
place_labels <- function(pos, count, half_width, marker_half, gap, dy) {
  y <- numeric(length(pos))
  boxes <- list()
  for (i in order(pos)) {
    x1 <- pos[i] - half_width[i]
    x2 <- pos[i] + half_width[i]
    cand <- count[i] + gap
    # clear any taller lollipop the label would otherwise sit on top of
    over <- which(pos + marker_half > x1 & pos - marker_half < x2)
    if (length(over)) cand <- max(cand, max(count[over]) + gap)
    repeat {
      hit <- FALSE
      for (b in boxes) {
        if (x1 < b$x2 && x2 > b$x1 && abs(cand - b$y) < dy * 0.95) {
          cand <- b$y + dy
          hit <- TRUE
          break
        }
      }
      if (!hit) break
    }
    boxes[[length(boxes) + 1L]] <- list(x1 = x1, x2 = x2, y = cand)
    y[i] <- cand
  }
  y
}

# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #

LABEL_SIZE <- 2.2   # ggplot size (mm) for mutation labels
GENE_SIZE  <- 2.5   # ggplot size (mm) for gene names

build_figure <- function(lps_type, chrom, locus_len, genes, variants, gcols, cfg) {
  xmax      <- max(locus_len, max(genes$end))
  max_count <- max(variants$count)
  n_var     <- nrow(variants)

  fig_w <- if (!is.null(cfg$width)) cfg$width else min(16, max(7, xmax / 1200 + n_var * 0.16))
  # rough panel width: figure minus the y axis furniture and right margin
  bp_per_inch <- xmax / max(2, fig_w - 1.1)

  # gene track lives below the count baseline, sized in count units
  band    <- max_count * 0.18
  arrow_h <- band / 2
  track_y <- -(band / 2) - max_count * 0.05

  genes$head <- pmin(0.35 * (genes$end - genes$start), 0.022 * xmax)
  genes$fill <- unname(gcols[genes$name])
  poly <- do.call(rbind, lapply(seq_len(nrow(genes)), function(i) {
    g <- genes[i, ]
    cbind(arrow_poly(g$start, g$end, track_y, arrow_h, g$strand, g$head),
          gene = g$name, id = i)
  }))
  poly$gene <- factor(poly$gene, levels = genes$name)

  # a gene name goes inside its arrow only if it fits in the straight body
  name_w <- measure_widths(genes$name, GENE_SIZE * ggplot2::.pt, font = 4L, family = cfg$font)
  genes$fits <- (name_w * bp_per_inch) < 0.9 * (genes$end - genes$start - genes$head)
  hl <- pmin(genes$head, 0.5 * (genes$end - genes$start))
  genes$mid  <- ifelse(genes$strand == "-",
                       (genes$start + hl + genes$end) / 2,
                       (genes$start + genes$end - hl) / 2)
  genes$text_colour <- contrast_text(genes$fill)

  inside  <- genes[genes$fits, , drop = FALSE]
  outside <- genes[!genes$fits, , drop = FALSE]

  # mutation labels, sitting above their own marker with leader lines
  lab_w <- measure_widths(variants$label, LABEL_SIZE * ggplot2::.pt, family = cfg$font)
  gap   <- max_count * 0.10
  dy    <- max_count * 0.15
  variants$y_lab <- place_labels(
    variants$pos, variants$count,
    half_width  = lab_w * bp_per_inch / 2 + xmax * 0.004,
    marker_half = xmax * 0.004,
    gap = gap, dy = dy)
  variants$fill <- factor(ifelse(variants$gene %in% genes$name, variants$gene, NA),
                          levels = genes$name)
  n_tiers <- max(1L, as.integer(ceiling((max(variants$y_lab) - max_count) / dy)))

  # height follows the count range and the label stack, so a run with only a
  # handful of genomes does not get a mostly empty panel
  fig_h <- if (!is.null(cfg$height)) cfg$height else
    2.3 + 0.9 * min(1, max_count / 8) + n_tiers * 0.22

  y_top     <- max(variants$y_lab) + max_count * 0.13
  arrow_bot <- track_y - arrow_h
  x_left    <- -0.02 * xmax  # panel edge, matching the x scale expansion below

  # Room under the track for gene names that did not fit inside their arrow.
  # Text height is physical (mm) while the y axis is in genome counts, so
  # reserve the space in inches and solve for the matching data range.
  if (nrow(outside)) {
    panel_in <- max(1.2, fig_h - 1.5)   # figure minus title, subtitle and x axis
    k <- min(0.35, (GENE_SIZE / 25.4 * 1.9) / panel_in)
    y_bottom <- arrow_bot - k * (y_top - arrow_bot) / (1 - k)
  } else {
    y_bottom <- arrow_bot - band * 0.35
  }
  gene_lab_y <- arrow_bot - 0.18 * (arrow_bot - y_bottom)

  step    <- max(1, ceiling(max_count / 5))
  ybreaks <- seq(0, max_count, by = step)

  leg_genes <- switch(cfg$legend,
                      all        = genes$name,
                      unlabelled = outside$name,
                      character(0))
  leg_labels <- sprintf("%s (%s-%s)",
                        genes$name[match(leg_genes, genes$name)],
                        scales::comma(genes$start[match(leg_genes, genes$name)]),
                        scales::comma(genes$end[match(leg_genes, genes$name)]))

  n_genomes <- attr(variants, "n_genomes")

  p <- ggplot() +
    # y axis spine, drawn by hand so it stops at the data instead of running
    # the full height of the label stack
    annotate("segment", x = x_left, xend = x_left, y = 0, yend = max_count,
             colour = "black", linewidth = 0.3) +
    # locus backbone, then arrows, then a single outline over the top
    annotate("segment", x = 1, xend = locus_len, y = track_y, yend = track_y,
             colour = "grey35", linewidth = 0.4) +
    geom_polygon(data = poly, aes(x = x, y = y, group = id, fill = gene), colour = NA) +
    geom_polygon(data = poly, aes(x = x, y = y, group = id),
                 fill = NA, colour = "grey20", linewidth = 0.3) +
    annotate("segment", x = x_left, xend = 1.02 * xmax, y = 0, yend = 0,
             colour = "grey60", linewidth = 0.4) +
    # lollipops
    geom_segment(data = variants, aes(x = pos, xend = pos, y = 0, yend = count),
                 colour = "grey75", linewidth = 0.4) +
    geom_segment(data = variants,
                 aes(x = pos, xend = pos, y = count, yend = y_lab - gap * 0.35),
                 colour = "grey85", linewidth = 0.3) +
    geom_point(data = variants, aes(x = pos, y = count, fill = fill),
               shape = 21, size = 2.6, colour = "white", stroke = 0.4,
               na.rm = TRUE) +
    geom_text(data = variants, aes(x = pos, y = y_lab, label = label),
              size = LABEL_SIZE, vjust = 0, colour = "grey10", family = cfg$font)

  if (nrow(inside)) {
    p <- p + geom_text(data = inside,
                       aes(x = mid, y = track_y, label = name, colour = text_colour),
                       size = GENE_SIZE, fontface = "bold.italic", family = cfg$font,
                       show.legend = FALSE)
  }
  if (nrow(outside)) {
    p <- p + geom_text(data = outside,
                       aes(x = mid, y = gene_lab_y, label = name),
                       size = GENE_SIZE, fontface = "bold.italic", vjust = 1,
                       colour = "grey20", family = cfg$font)
  }

  p <- p +
    scale_colour_identity() +
    scale_fill_manual(values = gcols, breaks = leg_genes, labels = leg_labels,
                      na.value = INTERGENIC_COLOUR, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = ybreaks, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = c(0, xmax), ylim = c(y_bottom, y_top), clip = "off") +
    labs(
      title    = sprintf("%s observed mutations", lps_type),
      subtitle = sprintf("%s, %s bp  \u00b7  %s genome%s with a called mutation  \u00b7  %s distinct mutation%s",
                         chrom, scales::comma(locus_len), n_genomes,
                         if (n_genomes == 1) "" else "s", n_var,
                         if (n_var == 1) "" else "s"),
      x = "position (bp)", y = "genomes with mutation") +
    theme_classic(base_size = 9, base_family = cfg$font) +
    theme(
      plot.title      = element_text(size = 11),
      plot.subtitle   = element_text(size = 7.5, colour = "grey35"),
      axis.title      = element_text(size = 8.5),
      axis.text       = element_text(size = 8),
      axis.line.y     = element_blank(),
      axis.line.x     = element_line(linewidth = 0.3),
      legend.position = if (length(leg_genes)) "bottom" else "none",
      legend.text     = element_text(size = 7.5),
      legend.key.size = unit(8, "pt"),
      legend.margin   = margin(0, 0, 0, 0),
      plot.margin     = margin(6, 12, 6, 6)
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  list(plot = p, width = fig_w, height = fig_h)
}

save_figure <- function(p, path, width, height, fmt, font) {
  dev <- switch(fmt,
    pdf = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
    svg = svglite::svglite,
    png = NULL)
  save <- function(plot) {
    if (is.null(dev)) {
      ggsave(path, plot, width = width, height = height, dpi = 600)
    } else {
      ggsave(path, plot, width = width, height = height, device = dev)
    }
  }
  tryCatch(save(p), error = function(e) {
    # most commonly an unavailable font family
    warning("could not render ", path, " with font '", font, "' (", conditionMessage(e),
            "); retrying with the device default", call. = FALSE)
    save(p + theme(text = element_text(family = "")))
  })
}

# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

main <- function() {
  db_dir <- opt("lps-db-dir")
  if (is.null(db_dir)) stop("--lps-db-dir is required\n", USAGE, call. = FALSE)
  if (!dir.exists(db_dir)) stop("no such directory: ", db_dir, call. = FALSE)

  reports <- if (!is.null(opt("subtype-report"))) {
    split_csv(opt("subtype-report"))
  } else {
    discover_reports(split_csv(opt("report-dir", ".")))
  }
  for (r in reports) if (!file.exists(r)) stop("no such file: ", r, call. = FALSE)
  message("reading: ", paste(reports, collapse = ", "))

  cfg <- list(
    outdir  = opt("outdir", "lps_mutation_figures"),
    prefix  = opt("prefix", "lps_mutations_"),
    formats = split_csv(opt("formats", "pdf,svg")),
    palette = opt("palette", "db"),
    legend  = opt("legend", "none"),
    font    = opt("font", "Helvetica"),
    width   = if (!is.null(opt("width"))) as.numeric(opt("width")) else NULL,
    height  = if (!is.null(opt("height"))) as.numeric(opt("height")) else NULL
  )
  if (!cfg$palette %in% c("db", "okabe-ito")) stop("--palette must be db or okabe-ito", call. = FALSE)
  if (!cfg$legend %in% c("all", "unlabelled", "none")) {
    stop("--legend must be all, unlabelled or none", call. = FALSE)
  }
  bad_fmt <- setdiff(cfg$formats, c("pdf", "svg", "png"))
  if (length(bad_fmt)) stop("unsupported format(s): ", paste(bad_fmt, collapse = ", "), call. = FALSE)

  min_genomes <- as.integer(opt("min-genomes", "1"))

  df <- read_subtype_reports(reports)

  # sample selection, applied before aggregation so counts match what is drawn
  keep <- split_csv(opt("samples"))
  if (!is.null(opt("samples-file"))) {
    keep <- c(keep, trimws(readLines(opt("samples-file"), warn = FALSE)))
  }
  keep <- unique(keep[nzchar(keep)])
  if (length(keep)) {
    unknown <- setdiff(keep, df$SAMPLE)
    if (length(unknown)) {
      warning("sample(s) not in the report: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    df <- df[df$SAMPLE %in% keep, , drop = FALSE]
  }
  drop <- split_csv(opt("exclude-samples"))
  if (length(drop)) df <- df[!df$SAMPLE %in% drop, , drop = FALSE]

  df <- df[!is_blank(df$POS) & !is.na(suppressWarnings(as.integer(df$POS))), , drop = FALSE]
  df <- df[!is_blank(df$TYPE) & df$TYPE != "untypeable", , drop = FALSE]
  if (!nrow(df)) stop("no mutation rows left after filtering", call. = FALSE)

  ref_map     <- read_reference_map(db_dir)
  colour_file <- opt("gene-colors", file.path(db_dir, "gene_colors.tsv"))
  gene_colors <- read_gene_colors(colour_file)

  want_types <- split_csv(opt("types"))
  types <- sort(unique(df$TYPE))
  if (length(want_types)) {
    unknown <- setdiff(want_types, types)
    if (length(unknown)) warning("no mutations for type(s): ", paste(unknown, collapse = ", "),
                                 call. = FALSE)
    types <- intersect(types, want_types)
  }
  if (!length(types)) stop("no LPS types to plot", call. = FALSE)

  dir.create(cfg$outdir, showWarnings = FALSE, recursive = TRUE)

  for (lps_type in types) {
    sub <- df[df$TYPE == lps_type, , drop = FALSE]
    gb_name <- ref_map[[lps_type]]
    if (is.null(gb_name)) { warning("no reference for ", lps_type, "; skipped", call. = FALSE); next }
    gb_path <- file.path(db_dir, gb_name)
    if (!file.exists(gb_path)) { warning("missing ", gb_path, "; skipped", call. = FALSE); next }

    ref <- parse_genbank_genes(gb_path)
    if (is.null(ref$genes)) { warning("no genes in ", gb_path, "; skipped", call. = FALSE); next }

    variants <- sub %>%
      mutate(pos = as.integer(POS), label = short_vartype(VARTYPE)) %>%
      group_by(CHROM, pos, REF, ALT) %>%
      summarise(gene = dplyr::first(GENE), label = dplyr::first(label),
                count = dplyr::n_distinct(SAMPLE), .groups = "drop") %>%
      filter(count >= min_genomes) %>%
      arrange(pos) %>%
      as.data.frame()
    if (!nrow(variants)) { warning("no mutations for ", lps_type, " after --min-genomes; skipped",
                                   call. = FALSE); next }
    variants$gene[is_blank(variants$gene)] <- NA_character_
    attr(variants, "n_genomes") <- length(unique(sub$SAMPLE))

    gcols <- color_for_genes(ref$genes$name, gene_colors, cfg$palette)
    chrom <- sub$CHROM[1]

    fig <- build_figure(lps_type, chrom, ref$length, ref$genes, variants, gcols, cfg)

    stem <- file.path(cfg$outdir, paste0(cfg$prefix, lps_type))
    for (fmt in cfg$formats) {
      out <- paste0(stem, ".", fmt)
      save_figure(fig$plot, out, fig$width, fig$height, fmt, cfg$font)
      message("wrote ", out)
    }
    if (isTRUE(opt("write-table"))) {
      out <- paste0(stem, ".tsv")
      utils::write.table(
        variants[c("CHROM", "pos", "REF", "ALT", "gene", "label", "count")],
        out, sep = "\t", quote = FALSE, row.names = FALSE)
      message("wrote ", out)
    }
  }
}

main()
