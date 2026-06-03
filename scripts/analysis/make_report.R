#!/usr/bin/env Rscript

# Generate Results-section plots and pivot tables from merged summaries.
#
# R/ggplot2 port of scripts/analysis/make_plots.py (avoids the matplotlib
# dependency). Same CLI flags, same output file names (TSVs + PNGs + README.md),
# so it is a drop-in replacement for the analysis_report rule.
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
  df$case <- paste0(tolower(df$mapper), ":", df[[label_col]])
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

CASE_PALETTE <- function(n) {
  # viridis discrete scale ships with ggplot2 and handles any number of levels.
  scale_fill_viridis_d(option = "D", begin = 0.05, end = 0.95)
}

grouped_bar <- function(df, index, columns, values, title, ylabel, path,
                        ylim = NULL, stacked = FALSE) {
  needed <- c(index, columns, values)
  d <- df[stats::complete.cases(df[, c(index, columns), drop = FALSE]), needed, drop = FALSE]
  if (nrow(d) == 0) return(invisible(NULL))
  d[[values]] <- as.numeric(d[[values]])
  agg <- aggregate(d[, values, drop = FALSE],
                   by = d[, c(index, columns), drop = FALSE],
                   FUN = function(x) mean(x, na.rm = TRUE))
  names(agg)[ncol(agg)] <- values
  if (all(is.na(agg[[values]]))) return(invisible(NULL))

  agg[[index]] <- factor(agg[[index]], levels = sort(unique(agg[[index]])))
  agg[[columns]] <- factor(agg[[columns]], levels = sort(unique(agg[[columns]])))

  position <- if (stacked) "stack" else position_dodge2(preserve = "single")
  plot <- ggplot(agg, aes(x = .data[[index]], y = .data[[values]], fill = .data[[columns]])) +
    geom_col(position = position, colour = "white", linewidth = 0.15, width = 0.85) +
    CASE_PALETTE() +
    labs(title = title, x = NULL, y = ylabel, fill = columns) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  if (!is.null(ylim)) plot <- plot + coord_cartesian(ylim = ylim)

  n_groups <- length(levels(agg[[index]]))
  width_in <- max(6, 0.35 * n_groups)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = width_in, height = 5, dpi = 150, bg = "white", limitsize = FALSE)
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
  df <- aggregate_precision(eligible, c("dataset", "mapper", "run", "as_event_type"))
  write_tsv(df, file.path(tables, "precision_by_as_event.tsv"))
  for (dataset in unique(df$dataset)) {
    subset <- with_label(df[df$dataset == dataset, , drop = FALSE])
    grouped_bar(subset, index = "case", columns = "as_event_type", values = "precision",
                title = sprintf("%s — precision by AS event", dataset),
                ylabel = "precision",
                path = file.path(plots, sprintf("precision_by_as_event_%s.png", dataset)),
                ylim = c(0, 1))
  }
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
    c("overall_alignment_rate", "HISAT2 overall alignment rate (%)")
  )
  available <- Filter(function(m) m[1] %in% names(qc), metric_cols)
  if (length(available) == 0) return(invisible(NULL))
  cols <- vapply(available, `[`, character(1), 1)
  summary <- group_means(qc, c("dataset", "mapper", "run"), cols)
  write_tsv(summary, file.path(tables, "mapper_qc.tsv"))
  for (metric in available) {
    col <- metric[1]; label <- metric[2]
    if (all(is.na(summary[[col]]))) next
    for (dataset in unique(summary$dataset)) {
      subset <- with_label(summary[summary$dataset == dataset, , drop = FALSE])
      if (all(is.na(subset[[col]]))) next
      grouped_bar(subset, index = "case", columns = "mapper", values = col,
                  title = sprintf("%s — %s (from mapper log)", dataset, label),
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
                path = file.path(plots, sprintf("rmats_truth_events_%s.png", dataset)))
    grouped_bar(subset, index = "case", columns = "event_type", values = "significant_events",
                title = sprintf("%s — rMATS significant events", dataset),
                ylabel = "# significant events",
                path = file.path(plots, sprintf("rmats_significant_events_%s.png", dataset)))
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
                stacked = TRUE)
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
