#!/usr/bin/env Rscript

# Generate Results-section plots and pivot tables from merged summaries.
# This is the R counterpart to make_plots.py, 
# using ggplot2 for plotting and base R for data manipulation.
# The two scripts produce the same outputs, but
#
# Inputs are the merged TSVs produced by the Snakemake pipeline:
#   --ground-summary       ground-truth per-sample summary
#   --stratified-summary   stratified ground-truth summary
#                          (read_type x transcript_type x as_event_type x error_type)
#   --standard-summary     standard mapping (BAM-derived QC) summary
#   --mapper-qc-summary    mapper-log QC summary (optional)
#   --rmats-summary        rMATS per-case summary (optional)
#   --outdir               output directory; plots/ and tables/ are created under it

suppressPackageStartupMessages(library(ggplot2))

AS_DATASETS <- c("sim_mouse_dataset2", "sim_mouse_dataset3")
ERROR_TYPES <- c("incorrect_position", "incorrect_position_within_region", "incorrect_junction")

# ----------------------------------------------------------------------------
# CLI parsing (base R, so ggplot2 is the only external dependency)
# ----------------------------------------------------------------------------

parse_args <- function(raw) {
  defaults <- list(
    `ground-summary` = NA_character_,
    `stratified-summary` = NA_character_,
    `standard-summary` = NA_character_,
    `mapper-qc-summary` = "",
    `rmats-summary` = "",
    outdir = NA_character_
  )
  i <- 1
  while (i <= length(raw)) {
    key <- raw[i]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    name <- sub("^--", "", key)
    if (!name %in% names(defaults)) stop(sprintf("Unknown flag: %s", key))
    if (i + 1 > length(raw)) stop(sprintf("Missing value for %s", key))
    defaults[[name]] <- raw[i + 1]
    i <- i + 2
  }
  required <- c("ground-summary", "stratified-summary", "standard-summary", "outdir")
  missing <- required[vapply(defaults[required], function(v) is.na(v) || v == "", logical(1))]
  if (length(missing)) stop(sprintf("Missing required flags: --%s", paste(missing, collapse = ", --")))
  defaults
}

# ----------------------------------------------------------------------------
# IO helpers
# ----------------------------------------------------------------------------

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

read_optional_tsv <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) return(NULL)
  df <- tryCatch(read_tsv(path), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

# ----------------------------------------------------------------------------
# Aggregation helpers (mirror the pandas groupby/sum/mean semantics)
# ----------------------------------------------------------------------------

with_label <- function(df, label_col = "run") {
  df$case <- df[[label_col]]
  df
}

# Sum count columns over the group, then derive precision = correct / (correct + incorrect).
aggregate_precision <- function(strat, group_cols) {
  count_cols <- c("n_reads", "n_correct", "n_incorrect", "n_unmapped")
  keep <- stats::complete.cases(strat[, group_cols, drop = FALSE])
  d <- strat[keep, , drop = FALSE]
  if (nrow(d) == 0) return(d[, c(group_cols, count_cols, "precision")][0, ])
  for (col in count_cols) d[[col]] <- as.numeric(d[[col]])
  agg <- aggregate(d[, count_cols, drop = FALSE], by = d[, group_cols, drop = FALSE],
                   FUN = function(x) sum(x, na.rm = TRUE))
  denom <- agg$n_correct + agg$n_incorrect
  agg$precision <- ifelse(denom > 0, agg$n_correct / denom, NA_real_)
  agg
}

# Mean of value_cols within each group (NaN if a group is entirely NA).
group_means <- function(df, group_cols, value_cols) {
  for (col in value_cols) df[[col]] <- as.numeric(df[[col]])
  aggregate(df[, value_cols, drop = FALSE], by = df[, group_cols, drop = FALSE],
            FUN = function(x) mean(x, na.rm = TRUE))
}

# Long-format reshape equivalent to DataFrame.melt.
melt_cols <- function(df, id_vars, value_vars, var_name, value_name) {
  parts <- lapply(value_vars, function(v) {
    sub <- df[, id_vars, drop = FALSE]
    sub[[var_name]] <- v
    sub[[value_name]] <- as.numeric(df[[v]])
    sub
  })
  do.call(rbind, parts)
}

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

PALETTE_COLORS <- c("#082848",  
                    "#5c140c",  
                    "#2a6a95",  
                    "#a43a06", 
                    "#052739", 
                    "#663205",
                    "#08306B")

category_palette <- function(n) {
  if (n <= length(PALETTE_COLORS)) return(PALETTE_COLORS[seq_len(n)])
  grDevices::colorRampPalette(PALETTE_COLORS)(n)
}

LEGEND_TITLES <- c(
  transcript_type = "Transcript structure",
  read_type = "Read type",
  as_event_type = "AS event type",
  as_event = "AS event type",
  error_type = "Error type",
  event_type = "AS event type",
  event_truth_class = "Simulation truth class",
  suppa2_match = "SUPPA2 event-type support",
  direction_match = "Direction vs. simulated truth",
  mapper = "Mapper",
  dataset = "Dataset"
)

legend_title <- function(col) {
  if (col %in% names(LEGEND_TITLES)) return(unname(LEGEND_TITLES[col]))
  tools::toTitleCase(gsub("_", " ", col))
}

# Tidy category labels shown in the legend: drop the "events_" prefix and turn
# underscores into spaces (e.g. "events_isoform_switch" -> "isoform switch").
pretty_levels <- function(x) gsub("_", " ", sub("^events_", "", x))

# Human-readable names for SUPPA2 AS-event types. Multi-type annotations (any
# value with a ";") collapse to a single "Multiple event types" group.
AS_EVENT_LABELS <- c(
  SE = "Skipped exon (SE)",
  MXE = "Mutually exclusive exons (MXE)",
  A5SS = "Alt. 5' splice site (A5SS)",
  A3SS = "Alt. 3' splice site (A3SS)",
  RI = "Retained intron (RI)",
  AF = "Alt. first exon (AF)",
  AL = "Alt. last exon (AL)",
  single_transcript = "Single transcript",
  same_intron_chain = "Same intron chain",
  complex_or_other = "Complex / other",
  A5SS_A3SS_or_MXE = "Complex / other"
)

as_event_label <- function(x) {
  out <- unname(AS_EVENT_LABELS[x])
  out[grepl(";", x)] <- "Multiple event types"
  out[is.na(out)] <- "Complex / other"
  out
}

# falls back to the device's default serif if the exact face is not installed.
THESIS_FONT <- "Times New Roman"
MIN_BAR_SLOT_IN <- 0.35  

grouped_bar <- function(df, index, columns, values, title, ylabel, path,
                        ylim = NULL, stacked = FALSE, width_scale = NULL,
                        facet = NULL, facet_ncol = 1, panel_height = 5.2) {
  group_keys <- c(index, columns, facet)  # facet is NULL when not faceting
  needed <- c(group_keys, values)
  d <- df[stats::complete.cases(df[, group_keys, drop = FALSE]), needed, drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  d[[values]] <- as.numeric(d[[values]])
  # Drop rows with no value (e.g. STAR-only / HISAT2-only QC metrics) so the
  # plot shows only populated bars instead of empty gaps.
  d <- d[is.finite(d[[values]]), , drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  agg <- aggregate(d[, values, drop = FALSE],
                   by = d[, group_keys, drop = FALSE],
                   FUN = function(x) mean(x, na.rm = TRUE))
  names(agg)[ncol(agg)] <- values
  if (all(is.na(agg[[values]]))) return(invisible(NULL))

  agg[[index]] <- factor(agg[[index]], levels = sort(unique(agg[[index]])))
  agg[[columns]] <- factor(agg[[columns]], levels = sort(unique(agg[[columns]])))
  n_fill <- nlevels(agg[[columns]])

  caption_text <- gsub("sim_mouse_dataset([0-9]+)", "Dataset \\1", title)

  position <- if (stacked) "stack" else position_dodge2(preserve = "single")
  plot <- ggplot(agg, aes(x = .data[[index]], y = .data[[values]], fill = .data[[columns]])) +
    geom_col(position = position, colour = "white", linewidth = 0.15, width = 0.85) +
    scale_fill_manual(values = category_palette(n_fill), labels = pretty_levels) +
    labs(caption = caption_text, x = NULL, y = ylabel, fill = legend_title(columns)) +
    theme_minimal(base_size = 12, base_family = THESIS_FONT) +
    theme(
      text = element_text(family = THESIS_FONT),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      plot.caption = element_text(hjust = 0.5, size = 12, margin = margin(t = 10)),
      plot.caption.position = "plot",
      legend.position = "right"
    )
  if (!is.null(ylim)) plot <- plot + coord_cartesian(ylim = ylim)

  if (is.null(facet)) {
    n_groups <- nlevels(agg[[index]])
    facet_rows <- 1
    facet_cols <- 1
  } else {
    # Split into panels (e.g. by mapper, or small multiples by event type).
    # mapper panels read better with upper-cased strip labels; leave already
    # human-readable facets (e.g. AS-event names) as-is.
    facet_labeller <- if (facet == "mapper") as_labeller(toupper) else "label_value"
    plot <- plot + facet_wrap(vars(.data[[facet]]), ncol = facet_ncol,
                              scales = "free_x", labeller = facet_labeller)
    per_facet <- tapply(as.character(agg[[index]]), agg[[facet]],
                        function(v) length(unique(v)))
    n_groups <- max(per_facet)
    n_facets <- length(per_facet)
    facet_cols <- min(facet_ncol, n_facets)
    facet_rows <- ceiling(n_facets / facet_ncol)
  }
  if (is.null(width_scale)) {
    # Default: size to the number of bars drawn in the busiest panel so each is
    # >= one MIN_BAR_SLOT_IN slot. Stacked bars are one-per-group; dodged groups
    # hold up to one bar per fill level present in the busiest group.
    if (stacked) {
      n_bars <- n_groups
    } else {
      slot_key <- if (is.null(facet)) as.character(agg[[index]]) else
        paste(agg[[facet]], agg[[index]])
      n_bars <- n_groups * max(table(slot_key))
    }
    panel_width <- MIN_BAR_SLOT_IN * n_bars
  } else {
    panel_width <- width_scale * n_groups  # explicit per-panel width override
  }
  width_in <- max(6, panel_width * facet_cols)
  height_in <- panel_height * facet_rows
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = width_in, height = height_in, dpi = 150,
         bg = "white", limitsize = FALSE)
}

# ----------------------------------------------------------------------------
# Report sections (one-to-one with make_plots.py)
# ----------------------------------------------------------------------------

precision_by_structure <- function(strat, plots, tables) {
  df <- aggregate_precision(strat, c("dataset", "mapper", "run", "transcript_type"))
  write_tsv(df, file.path(tables, "precision_by_structure.tsv"))
  for (dataset in unique(df$dataset)) {
    subset <- with_label(df[df$dataset == dataset, , drop = FALSE])
    grouped_bar(subset, index = "case", columns = "transcript_type", values = "precision",
                title = sprintf("%s — precision by transcript structure", dataset),
                ylabel = "precision",
                path = file.path(plots, sprintf("precision_by_structure_%s.png", dataset)),
                ylim = c(0, 1))
  }
}

precision_by_read_type <- function(strat, plots, tables) {
  df <- aggregate_precision(strat, c("dataset", "mapper", "run", "read_type"))
  write_tsv(df, file.path(tables, "precision_by_read_type.tsv"))
  for (dataset in unique(df$dataset)) {
    subset <- with_label(df[df$dataset == dataset, , drop = FALSE])
    grouped_bar(subset, index = "case", columns = "read_type", values = "precision",
                title = sprintf("%s — precision by read type", dataset),
                ylabel = "precision",
                path = file.path(plots, sprintf("precision_by_read_type_%s.png", dataset)),
                ylim = c(0, 1))
  }
}

error_breakdown <- function(strat, plots, tables) {
  err <- strat[strat$error_type %in% ERROR_TYPES, , drop = FALSE]
  if (nrow(err) == 0) return(invisible(NULL))
  err$n_incorrect <- as.numeric(err$n_incorrect)
  agg <- aggregate(n_incorrect ~ dataset + mapper + run + error_type, err, FUN = function(x) sum(x, na.rm = TRUE))
  totals <- aggregate(n_incorrect ~ dataset + mapper + run, agg, FUN = function(x) sum(x, na.rm = TRUE))
  names(totals)[names(totals) == "n_incorrect"] <- "total_incorrect"
  agg <- merge(agg, totals, by = c("dataset", "mapper", "run"))
  agg$pct_of_incorrect <- ifelse(agg$total_incorrect > 0,
                                 100 * agg$n_incorrect / agg$total_incorrect, NA_real_)
  write_tsv(agg, file.path(tables, "error_breakdown.tsv"))
  for (dataset in unique(agg$dataset)) {
    subset <- with_label(agg[agg$dataset == dataset, , drop = FALSE])
    grouped_bar(subset, index = "case", columns = "error_type", values = "pct_of_incorrect",
                title = sprintf("%s — error-type breakdown", dataset),
                ylabel = "% of incorrect alignments",
                path = file.path(plots, sprintf("error_breakdown_%s.png", dataset)),
                stacked = TRUE)
  }
}

precision_by_as_event <- function(strat, plots, tables) {
  eligible <- strat[strat$dataset %in% AS_DATASETS, , drop = FALSE]
  if (nrow(eligible) == 0) return(invisible(NULL))
  eligible$as_event <- as_event_label(eligible$as_event_type)
  df <- aggregate_precision(eligible, c("dataset", "mapper", "as_event"))
  write_tsv(df, file.path(tables, "precision_by_as_event.tsv"))
  lo <- max(0, floor(min(df$precision, na.rm = TRUE) / 0.05) * 0.05 - 0.05)
  for (dataset in unique(df$dataset)) {
    subset <- df[df$dataset == dataset, , drop = FALSE]
    grouped_bar(subset, index = "as_event", columns = "mapper", values = "precision",
                title = sprintf("%s — precision by AS event", dataset),
                ylabel = "precision",
                path = file.path(plots, sprintf("precision_by_as_event_%s.png", dataset)),
                ylim = c(lo, 1))
  }
}

# Compare precision across all datasets on a shared axis. The dataset token
# (_d1/_d2/_d3) is stripped from the run so the same mapper-parameter setting
# (e.g. "hisat2_mm10_int100") lines up across datasets, with dataset as the fill.
precision_across_datasets <- function(strat, plots, tables) {
  d <- strat
  d$param_setting <- sub("_d[0-9]+", "", d$run)
  df <- aggregate_precision(d, c("param_setting", "dataset"))
  df$dataset <- sub("sim_mouse_dataset", "Dataset ", df$dataset)
  write_tsv(df, file.path(tables, "precision_across_datasets.tsv"))
  grouped_bar(df, index = "param_setting", columns = "dataset", values = "precision",
              title = "precision by parameter setting across datasets",
              ylabel = "precision",
              path = file.path(plots, "precision_across_datasets.png"),
              ylim = c(0, 1))
}

standard_qc <- function(standard, plots, tables) {
  df <- standard
  df$total_reads <- as.numeric(df$total_reads)
  df$unmapped_reads <- as.numeric(df$unmapped_reads)
  df$reads_with_multiple_alignments <- as.numeric(df$reads_with_multiple_alignments)
  df$pct_unmapped <- ifelse(df$total_reads > 0, 100 * df$unmapped_reads / df$total_reads, NA_real_)
  df$pct_multi <- ifelse(df$total_reads > 0, 100 * df$reads_with_multiple_alignments / df$total_reads, NA_real_)
  summary <- group_means(df, c("dataset", "mapper", "run"), c("pct_mapped", "pct_unmapped", "pct_multi"))
  write_tsv(summary, file.path(tables, "standard_qc.tsv"))
  metrics <- list(c("pct_mapped", "% mapped"),
                  c("pct_unmapped", "% unmapped"),
                  c("pct_multi", "% multi-aligned"))
  for (metric in metrics) {
    col <- metric[1]; label <- metric[2]
    for (dataset in unique(summary$dataset)) {
      subset <- with_label(summary[summary$dataset == dataset, , drop = FALSE])
      grouped_bar(subset, index = "case", columns = "mapper", values = col,
                  title = sprintf("%s — %s", dataset, label), ylabel = label,
                  path = file.path(plots, sprintf("standard_%s_%s.png", col, dataset)))
    }
  }
}

mapper_qc <- function(qc, plots, tables) {
  if (is.null(qc)) return(invisible(NULL))
  metric_cols <- list(
    c("pct_uniquely_mapped", "% uniquely mapped"),
    c("pct_multi_mapped", "% multi-mapped"),
    c("pct_unmapped", "% unmapped"),
    c("mismatch_rate_pct", "mismatch rate (%)"),
    c("overall_alignment_rate", "overall alignment rate (%)")
  )
  available <- Filter(function(m) m[1] %in% names(qc), metric_cols)
  if (length(available) == 0) return(invisible(NULL))
  cols <- vapply(available, `[`, character(1), 1)
  summary <- group_means(qc, c("dataset", "mapper", "run"), cols)
  write_tsv(summary, file.path(tables, "mapper_qc.tsv"))
  for (metric in available) {
    col <- metric[1]; label <- metric[2]
    if (all(is.na(summary[[col]]))) next
    # Some metrics are reported by only one mapper (mismatch rate = STAR,
    # overall alignment rate = HISAT2). Name that mapper in the title so the
    # single-colour plot is self-explanatory.
    mappers_present <- unique(summary$mapper[is.finite(summary[[col]])])
    only <- if (length(mappers_present) == 1) {
      sprintf(" (%s only)", toupper(mappers_present))
    } else {
      ""
    }
    for (dataset in unique(summary$dataset)) {
      subset <- with_label(summary[summary$dataset == dataset, , drop = FALSE])
      if (all(is.na(subset[[col]]))) next
      grouped_bar(subset, index = "case", columns = "mapper", values = col,
                  title = sprintf("%s — %s (from mapper log)%s", dataset, label, only),
                  ylabel = label,
                  path = file.path(plots, sprintf("mapper_qc_%s_%s.png", col, dataset)))
    }
  }
}

rmats_recovery <- function(rmats, plots, tables) {
  if (is.null(rmats)) return(invisible(NULL))
  df <- rmats
  df$significant_genes <- as.numeric(df$significant_genes)
  df$truth_supported_genes <- as.numeric(df$truth_supported_genes)
  df$significant_events <- as.numeric(df$significant_events)
  df$truth_supported_events <- as.numeric(df$truth_supported_events)
  df$gene_recovery_rate <- ifelse(df$significant_genes > 0,
                                  df$truth_supported_genes / df$significant_genes, NA_real_)
  df$event_recovery_rate <- ifelse(df$significant_events > 0,
                                   df$truth_supported_events / df$significant_events, NA_real_)
  write_tsv(df, file.path(tables, "rmats_truth_recovery.tsv"))
  for (dataset in unique(df$dataset)) {
    subset <- with_label(df[df$dataset == dataset, , drop = FALSE], "param_group")
    grouped_bar(subset, index = "case", columns = "event_type", values = "truth_supported_events",
                title = sprintf("%s — rMATS truth-supported events", dataset),
                ylabel = "# truth-supported events",
                path = file.path(plots, sprintf("rmats_truth_events_%s.png", dataset)),
                facet = "mapper")
    grouped_bar(subset, index = "case", columns = "event_type", values = "significant_events",
                title = sprintf("%s — rMATS significant events", dataset),
                ylabel = "# significant events",
                path = file.path(plots, sprintf("rmats_significant_events_%s.png", dataset)),
                facet = "mapper")
  }
}

# Shared body for the three rMATS class-composition stacked-bar sections.
rmats_class_composition <- function(rmats, plots, tables, class_cols, var_name,
                                    table_name, title_suffix, plot_prefix) {
  if (is.null(rmats)) return(invisible(NULL))
  if (!all(class_cols %in% names(rmats))) return(invisible(NULL))
  id_vars <- c("dataset", "mapper", "param_group", "event_type")
  long <- melt_cols(rmats, id_vars, class_cols, var_name, "n_events")
  write_tsv(long, file.path(tables, table_name))
  group_cols <- c("dataset", "mapper", "param_group", var_name)
  per_case <- aggregate(long[, "n_events", drop = FALSE], by = long[, group_cols, drop = FALSE],
                        FUN = function(x) sum(x, na.rm = TRUE))
  for (dataset in unique(per_case$dataset)) {
    subset <- per_case[per_case$dataset == dataset, , drop = FALSE]
    if (nrow(subset) == 0) next
    subset <- with_label(subset, "param_group")
    grouped_bar(subset, index = "case", columns = var_name, values = "n_events",
                title = sprintf("%s — %s", dataset, title_suffix),
                ylabel = "# significant events",
                path = file.path(plots, sprintf("%s_%s.png", plot_prefix, dataset)),
                stacked = TRUE, facet = "mapper")
  }
}

rmats_event_class_breakdown <- function(rmats, plots, tables) {
  rmats_class_composition(
    rmats, plots, tables,
    class_cols = c("events_isoform_switch", "events_single_shift",
                   "events_co_directional", "events_null_gene", "events_unknown_gene"),
    var_name = "event_truth_class",
    table_name = "rmats_event_class_breakdown.tsv",
    title_suffix = "rMATS significant events by truth class",
    plot_prefix = "rmats_event_class")
}

rmats_suppa2_type_match <- function(rmats, plots, tables) {
  rmats_class_composition(
    rmats, plots, tables,
    class_cols = c("events_suppa2_type_present", "events_suppa2_type_absent",
                   "events_suppa2_no_annotation"),
    var_name = "suppa2_match",
    table_name = "rmats_suppa2_type_match.tsv",
    title_suffix = "rMATS events vs SUPPA2 gene-level event types",
    plot_prefix = "rmats_suppa2_type_match")
}

rmats_direction_match <- function(rmats, plots, tables) {
  rmats_class_composition(
    rmats, plots, tables,
    class_cols = c("events_direction_match", "events_direction_mismatch",
                   "events_direction_indeterminate", "events_direction_na"),
    var_name = "direction_match",
    table_name = "rmats_direction_match.tsv",
    title_suffix = "rMATS direction match vs simulation truth",
    plot_prefix = "rmats_direction_match")
}

write_index <- function(outdir) {
  lines <- c("# Pipeline analysis outputs", "",
             "Generated by `scripts/analysis/make_report.R`.", "")
  for (sub in c("plots", "tables")) {
    section_dir <- file.path(outdir, sub)
    if (!dir.exists(section_dir)) next
    lines <- c(lines, sprintf("## %s", tools::toTitleCase(sub)), "")
    for (name in sort(list.files(section_dir))) {
      lines <- c(lines, sprintf("- `%s/%s`", sub, name))
    }
    lines <- c(lines, "")
  }
  writeLines(lines, file.path(outdir, "README.md"))
}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  outdir <- args$outdir
  plots <- file.path(outdir, "plots")
  tables <- file.path(outdir, "tables")
  dir.create(plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  strat <- read_tsv(args$`stratified-summary`)
  standard <- read_tsv(args$`standard-summary`)
  qc <- read_optional_tsv(args$`mapper-qc-summary`)
  rmats <- read_optional_tsv(args$`rmats-summary`)

  precision_by_structure(strat, plots, tables)
  precision_by_read_type(strat, plots, tables)
  precision_across_datasets(strat, plots, tables)
  error_breakdown(strat, plots, tables)
  precision_by_as_event(strat, plots, tables)
  standard_qc(standard, plots, tables)
  mapper_qc(qc, plots, tables)
  rmats_recovery(rmats, plots, tables)
  rmats_event_class_breakdown(rmats, plots, tables)
  rmats_suppa2_type_match(rmats, plots, tables)
  rmats_direction_match(rmats, plots, tables)

  write_index(outdir)
}

main()
