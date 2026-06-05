#!/usr/bin/env Rscript
# Generate Results-section plots and pivot tables from merged summaries.
#
# Inputs (merged TSVs from the Snakemake pipeline):
#   --ground-summary       ground-truth per-sample summary
#   --stratified-summary   stratified ground-truth summary
#   --standard-summary     standard mapping (BAM-derived QC) summary
#   --mapper-qc-summary    mapper-log QC summary (optional)
#   --rmats-summary        rMATS per-case summary (optional)
#   --outdir               output directory; plots/ and tables/ created under it

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

AS_DATASETS <- c("sim_mouse_dataset2", "sim_mouse_dataset3")
ERROR_TYPES <- c("incorrect_position", "incorrect_position_within_region", "incorrect_junction")

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

parse_args <- function(raw) {
  defaults <- list(
    `ground-summary`     = NA_character_,
    `stratified-summary` = NA_character_,
    `standard-summary`   = NA_character_,
    `mapper-qc-summary`  = "",
    `rmats-summary`      = "",
    outdir               = NA_character_
  )
  i <- 1
  while (i <= length(raw)) {
    key  <- raw[i]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    name <- sub("^--", "", key)
    if (!name %in% names(defaults)) stop(sprintf("Unknown flag: %s", key))
    if (i + 1 > length(raw)) stop(sprintf("Missing value for %s", key))
    defaults[[name]] <- raw[i + 1]
    i <- i + 2
  }
  required <- c("ground-summary", "stratified-summary", "standard-summary", "outdir")
  missing  <- required[vapply(defaults[required], function(v) is.na(v) || v == "", logical(1))]
  if (length(missing)) stop(sprintf("Missing required flags: --%s", paste(missing, collapse = ", --")))
  defaults
}

# ----------------------------------------------------------------------------
# IO
# ----------------------------------------------------------------------------

read_tsv <- function(path) {
  read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
             stringsAsFactors = FALSE, na.strings = c("NA", ""))
}

read_optional_tsv <- function(path) {
  if (is.null(path) || is.na(path) || path == "" || !file.exists(path)) return(NULL)
  tryCatch({ df <- read_tsv(path); if (nrow(df) == 0) NULL else df }, error = function(e) NULL)
}

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

# ----------------------------------------------------------------------------
# Aggregation helpers
# ----------------------------------------------------------------------------

aggregate_precision <- function(df, group_cols) {
  count_cols <- c("n_reads", "n_correct", "n_incorrect", "n_unmapped")
  df |>
    filter(if_all(all_of(group_cols), ~ !is.na(.x))) |>
    mutate(across(all_of(count_cols), as.numeric)) |>
    group_by(across(all_of(group_cols))) |>
    summarise(across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop") |>
    mutate(precision = if_else(n_correct + n_incorrect > 0,
                               n_correct / (n_correct + n_incorrect), NA_real_))
}

group_means <- function(df, group_cols, value_cols) {
  df |>
    mutate(across(all_of(value_cols), as.numeric)) |>
    group_by(across(all_of(group_cols))) |>
    summarise(across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

PALETTE_COLORS <- c("#08519C", "#70160c", "#1b374b", "#5a250b", "#3c5f71", "#9d591d", "#08306B")

category_palette <- function(n) {
  if (n <= length(PALETTE_COLORS)) PALETTE_COLORS[seq_len(n)]
  else grDevices::colorRampPalette(PALETTE_COLORS)(n)
}

LEGEND_TITLES <- c(
  transcript_type   = "Transcript structure",
  read_type         = "Read type",
  as_event_type     = "AS event type",
  as_event          = "AS event type",
  error_type        = "Error type",
  event_type        = "AS event type",
  event_truth_class = "Simulation truth class",
  suppa2_match      = "SUPPA2 event-type support",
  direction_match   = "Direction vs. simulated truth",
  mapper            = "Mapper",
  dataset           = "Dataset"
)

legend_title  <- function(col) {
  if (col %in% names(LEGEND_TITLES)) unname(LEGEND_TITLES[col])
  else tools::toTitleCase(gsub("_", " ", col))
}
pretty_levels <- function(x) gsub("_", " ", sub("^events_", "", x))

AS_EVENT_LABELS <- c(
  SE = "Skipped exon (SE)", MXE = "Mutually exclusive exons (MXE)",
  A5SS = "Alt. 5' splice site (A5SS)", A3SS = "Alt. 3' splice site (A3SS)",
  RI = "Retained intron (RI)", AF = "Alt. first exon (AF)", AL = "Alt. last exon (AL)",
  single_transcript = "Single transcript", same_intron_chain = "Same intron chain",
  complex_or_other = "Complex / other", A5SS_A3SS_or_MXE = "Complex / other"
)

as_event_label <- function(x) {
  out <- unname(AS_EVENT_LABELS[x])
  out[grepl(";", x)] <- "Multiple event types"
  out[is.na(out)]    <- "Complex / other"
  out
}

THESIS_FONT  <- "Times New Roman"
MIN_BAR_SLOT <- 0.35

grouped_bar <- function(df, index, columns, values, title, ylabel, path,
                        ylim = NULL, stacked = FALSE, width_scale = NULL,
                        facet = NULL, facet_ncol = 1, panel_height = 5.2) {
  group_keys <- c(index, columns, facet)
  d <- df |>
    filter(if_all(all_of(group_keys), ~ !is.na(.x))) |>
    mutate(across(all_of(values), as.numeric)) |>
    filter(is.finite(.data[[values]]))
  if (nrow(d) == 0) return(invisible(NULL))

  agg <- d |>
    group_by(across(all_of(group_keys))) |>
    summarise(across(all_of(values), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  if (all(is.na(agg[[values]]))) return(invisible(NULL))

  agg[[index]]   <- factor(agg[[index]],   levels = sort(unique(agg[[index]])))
  agg[[columns]] <- factor(agg[[columns]], levels = sort(unique(agg[[columns]])))
  n_fill         <- nlevels(agg[[columns]])
  caption_text   <- gsub("sim_mouse_dataset([0-9]+)", "Dataset \\1", title)
  position       <- if (stacked) "stack" else position_dodge2(preserve = "single")

  plot <- ggplot(agg, aes(x = .data[[index]], y = .data[[values]], fill = .data[[columns]])) +
    geom_col(position = position, colour = "white", linewidth = 0.15, width = 0.85) +
    scale_fill_manual(values = category_palette(n_fill), labels = pretty_levels) +
    labs(caption = caption_text, x = NULL, y = ylabel, fill = legend_title(columns)) +
    theme_minimal(base_size = 12, base_family = THESIS_FONT) +
    theme(
      text                  = element_text(family = THESIS_FONT),
      axis.text.x           = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
      axis.text.y           = element_text(size = 10),
      panel.grid.major.x    = element_blank(),
      plot.caption          = element_text(hjust = 0.5, size = 12, margin = margin(t = 10)),
      plot.caption.position = "plot",
      legend.position       = "right"
    )
  if (!is.null(ylim)) plot <- plot + coord_cartesian(ylim = ylim)

  if (is.null(facet)) {
    n_groups <- nlevels(agg[[index]]); facet_rows <- 1; facet_cols <- 1
  } else {
    labeller  <- if (facet == "mapper") as_labeller(toupper) else "label_value"
    plot      <- plot + facet_wrap(vars(.data[[facet]]), ncol = facet_ncol,
                                   scales = "free_x", labeller = labeller)
    per_facet <- tapply(as.character(agg[[index]]), agg[[facet]], function(v) length(unique(v)))
    n_groups   <- max(per_facet)
    facet_cols <- min(facet_ncol, length(per_facet))
    facet_rows <- ceiling(length(per_facet) / facet_ncol)
  }

  slot_key    <- if (is.null(facet) || stacked) as.character(agg[[index]]) else
    paste(agg[[facet]], agg[[index]])
  n_bars      <- if (stacked) n_groups else n_groups * max(table(slot_key))
  panel_width <- if (is.null(width_scale)) MIN_BAR_SLOT * n_bars else width_scale * n_groups

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(path, plot, width = max(6, panel_width * facet_cols),
         height = panel_height * facet_rows, dpi = 150, bg = "white", limitsize = FALSE)
}

# ----------------------------------------------------------------------------
# Report sections
# ----------------------------------------------------------------------------

precision_by_structure <- function(strat, plots, tables) {
  df <- aggregate_precision(strat, c("dataset", "mapper", "run", "transcript_type"))
  write_tsv(df, file.path(tables, "precision_by_structure.tsv"))
  for (ds in unique(df$dataset)) {
    grouped_bar(mutate(filter(df, dataset == ds), case = run),
      index = "case", columns = "transcript_type", values = "precision",
      title = sprintf("%s — precision by transcript structure", ds), ylabel = "precision",
      path  = file.path(plots, sprintf("precision_by_structure_%s.png", ds)), ylim = c(0, 1))
  }
}

precision_by_read_type <- function(strat, plots, tables) {
  df <- aggregate_precision(strat, c("dataset", "mapper", "run", "read_type"))
  write_tsv(df, file.path(tables, "precision_by_read_type.tsv"))
  for (ds in unique(df$dataset)) {
    grouped_bar(mutate(filter(df, dataset == ds), case = run),
      index = "case", columns = "read_type", values = "precision",
      title = sprintf("%s — precision by read type", ds), ylabel = "precision",
      path  = file.path(plots, sprintf("precision_by_read_type_%s.png", ds)), ylim = c(0, 1))
  }
}

precision_across_datasets <- function(strat, plots, tables) {
  df <- strat |>
    mutate(param_setting = sub("_d[0-9]+", "", run),
           dataset       = sub("sim_mouse_dataset", "Dataset ", dataset)) |>
    aggregate_precision(c("param_setting", "dataset"))
  write_tsv(df, file.path(tables, "precision_across_datasets.tsv"))
  grouped_bar(df, index = "param_setting", columns = "dataset", values = "precision",
    title = "precision by parameter setting across datasets", ylabel = "precision",
    path  = file.path(plots, "precision_across_datasets.png"), ylim = c(0, 1))
}

error_breakdown <- function(strat, plots, tables) {
  agg <- strat |>
    filter(error_type %in% ERROR_TYPES) |>
    mutate(n_incorrect = as.numeric(n_incorrect)) |>
    group_by(dataset, mapper, run, error_type) |>
    summarise(n_incorrect = sum(n_incorrect, na.rm = TRUE), .groups = "drop") |>
    group_by(dataset, mapper, run) |>
    mutate(pct_of_incorrect = 100 * n_incorrect / sum(n_incorrect)) |>
    ungroup()
  if (nrow(agg) == 0) return(invisible(NULL))
  write_tsv(agg, file.path(tables, "error_breakdown.tsv"))
  for (ds in unique(agg$dataset)) {
    grouped_bar(mutate(filter(agg, dataset == ds), case = run),
      index = "case", columns = "error_type", values = "pct_of_incorrect",
      title = sprintf("%s — error-type breakdown", ds), ylabel = "% of incorrect alignments",
      path  = file.path(plots, sprintf("error_breakdown_%s.png", ds)), stacked = TRUE)
  }
}

precision_by_as_event <- function(strat, plots, tables) {
  eligible <- strat |>
    filter(dataset %in% AS_DATASETS) |>
    mutate(as_event = as_event_label(as_event_type))
  if (nrow(eligible) == 0) return(invisible(NULL))
  df <- aggregate_precision(eligible, c("dataset", "mapper", "as_event"))
  write_tsv(df, file.path(tables, "precision_by_as_event.tsv"))
  lo <- max(0, floor(min(df$precision, na.rm = TRUE) / 0.05) * 0.05 - 0.05)
  for (ds in unique(df$dataset)) {
    grouped_bar(filter(df, dataset == ds),
      index = "as_event", columns = "mapper", values = "precision",
      title = sprintf("%s — precision by AS event", ds), ylabel = "precision",
      path  = file.path(plots, sprintf("precision_by_as_event_%s.png", ds)), ylim = c(lo, 1))
  }
}

standard_qc <- function(standard, plots, tables) {
  df <- standard |>
    mutate(across(c(total_reads, unmapped_reads, reads_with_multiple_alignments), as.numeric),
           pct_unmapped = if_else(total_reads > 0, 100 * unmapped_reads / total_reads, NA_real_),
           pct_multi    = if_else(total_reads > 0,
                                  100 * reads_with_multiple_alignments / total_reads, NA_real_))
  summary <- group_means(df, c("dataset", "mapper", "run"),
                         c("pct_mapped", "pct_unmapped", "pct_multi"))
  write_tsv(summary, file.path(tables, "standard_qc.tsv"))
  metrics <- list(c("pct_mapped", "% mapped"), c("pct_unmapped", "% unmapped"),
                  c("pct_multi", "% multi-aligned"))
  for (m in metrics) {
    for (ds in unique(summary$dataset)) {
      grouped_bar(mutate(filter(summary, dataset == ds), case = run),
        index = "case", columns = "mapper", values = m[1],
        title = sprintf("%s — %s", ds, m[2]), ylabel = m[2],
        path  = file.path(plots, sprintf("standard_%s_%s.png", m[1], ds)))
    }
  }
}

mapper_qc <- function(qc, plots, tables) {
  if (is.null(qc)) return(invisible(NULL))
  metrics <- list(
    c("pct_uniquely_mapped",    "% uniquely mapped"),
    c("pct_multi_mapped",       "% multi-mapped"),
    c("pct_unmapped",           "% unmapped"),
    c("mismatch_rate_pct",      "mismatch rate (%)"),
    c("overall_alignment_rate", "overall alignment rate (%)")
  )
  available <- Filter(function(m) m[1] %in% names(qc), metrics)
  if (length(available) == 0) return(invisible(NULL))
  cols    <- vapply(available, `[`, character(1), 1)
  summary <- group_means(qc, c("dataset", "mapper", "run"), cols)
  write_tsv(summary, file.path(tables, "mapper_qc.tsv"))
  for (m in available) {
    col <- m[1]; label <- m[2]
    if (all(is.na(summary[[col]]))) next
    mappers <- unique(summary$mapper[is.finite(summary[[col]])])
    only    <- if (length(mappers) == 1) sprintf(" (%s only)", toupper(mappers)) else ""
    for (ds in unique(summary$dataset)) {
      sub <- mutate(filter(summary, dataset == ds), case = run)
      if (all(is.na(sub[[col]]))) next
      grouped_bar(sub, index = "case", columns = "mapper", values = col,
        title  = sprintf("%s — %s (from mapper log)%s", ds, label, only), ylabel = label,
        path   = file.path(plots, sprintf("mapper_qc_%s_%s.png", col, ds)))
    }
  }
}

rmats_recovery <- function(rmats, plots, tables) {
  if (is.null(rmats)) return(invisible(NULL))
  df <- rmats |>
    mutate(across(c(significant_genes, truth_supported_genes,
                    significant_events, truth_supported_events), as.numeric),
           gene_recovery_rate  = if_else(significant_genes  > 0,
                                         truth_supported_genes  / significant_genes,  NA_real_),
           event_recovery_rate = if_else(significant_events > 0,
                                         truth_supported_events / significant_events, NA_real_))
  write_tsv(df, file.path(tables, "rmats_truth_recovery.tsv"))
  for (ds in unique(df$dataset)) {
    sub <- mutate(filter(df, dataset == ds), case = param_group)
    grouped_bar(sub, index = "case", columns = "event_type", values = "truth_supported_events",
      title = sprintf("%s — rMATS truth-supported events", ds), ylabel = "# truth-supported events",
      path  = file.path(plots, sprintf("rmats_truth_events_%s.png", ds)), facet = "mapper")
    grouped_bar(sub, index = "case", columns = "event_type", values = "significant_events",
      title = sprintf("%s — rMATS significant events", ds), ylabel = "# significant events",
      path  = file.path(plots, sprintf("rmats_significant_events_%s.png", ds)), facet = "mapper")
  }
}

rmats_class_composition <- function(rmats, plots, tables, class_cols, var_name,
                                    table_name, title_suffix, plot_prefix) {
  if (is.null(rmats) || !all(class_cols %in% names(rmats))) return(invisible(NULL))
  long <- rmats |>
    select(dataset, mapper, param_group, event_type, all_of(class_cols)) |>
    pivot_longer(all_of(class_cols), names_to = var_name, values_to = "n_events")
  write_tsv(long, file.path(tables, table_name))
  per_case <- long |>
    group_by(across(all_of(c("dataset", "mapper", "param_group", var_name)))) |>
    summarise(n_events = sum(as.numeric(n_events), na.rm = TRUE), .groups = "drop")
  for (ds in unique(per_case$dataset)) {
    grouped_bar(mutate(filter(per_case, dataset == ds), case = param_group),
      index = "case", columns = var_name, values = "n_events",
      title = sprintf("%s — %s", ds, title_suffix), ylabel = "# significant events",
      path  = file.path(plots, sprintf("%s_%s.png", plot_prefix, ds)),
      stacked = TRUE, facet = "mapper")
  }
}

rmats_event_class_breakdown <- function(rmats, plots, tables) {
  rmats_class_composition(rmats, plots, tables,
    class_cols   = c("events_isoform_switch", "events_single_shift",
                     "events_co_directional", "events_null_gene", "events_unknown_gene"),
    var_name     = "event_truth_class",
    table_name   = "rmats_event_class_breakdown.tsv",
    title_suffix = "rMATS significant events by truth class",
    plot_prefix  = "rmats_event_class")
}

rmats_suppa2_type_match <- function(rmats, plots, tables) {
  rmats_class_composition(rmats, plots, tables,
    class_cols   = c("events_suppa2_type_present", "events_suppa2_type_absent",
                     "events_suppa2_no_annotation"),
    var_name     = "suppa2_match",
    table_name   = "rmats_suppa2_type_match.tsv",
    title_suffix = "rMATS events vs SUPPA2 gene-level event types",
    plot_prefix  = "rmats_suppa2_type_match")
}

rmats_direction_match <- function(rmats, plots, tables) {
  rmats_class_composition(rmats, plots, tables,
    class_cols   = c("events_direction_match", "events_direction_mismatch",
                     "events_direction_indeterminate", "events_direction_na"),
    var_name     = "direction_match",
    table_name   = "rmats_direction_match.tsv",
    title_suffix = "rMATS direction match vs simulation truth",
    plot_prefix  = "rmats_direction_match")
}

write_index <- function(outdir) {
  lines <- c("# Pipeline analysis outputs", "",
             "Generated by `scripts/analysis/make_report.R`.", "")
  for (sub in c("plots", "tables")) {
    section_dir <- file.path(outdir, sub)
    if (!dir.exists(section_dir)) next
    lines <- c(lines, sprintf("## %s", tools::toTitleCase(sub)), "")
    for (name in sort(list.files(section_dir)))
      lines <- c(lines, sprintf("- `%s/%s`", sub, name))
    lines <- c(lines, "")
  }
  writeLines(lines, file.path(outdir, "README.md"))
}

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

main <- function() {
  args   <- parse_args(commandArgs(trailingOnly = TRUE))
  plots  <- file.path(args$outdir, "plots")
  tables <- file.path(args$outdir, "tables")
  dir.create(plots,  recursive = TRUE, showWarnings = FALSE)
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)

  strat    <- read_tsv(args$`stratified-summary`)
  standard <- read_tsv(args$`standard-summary`)
  qc       <- read_optional_tsv(args$`mapper-qc-summary`)
  rmats    <- read_optional_tsv(args$`rmats-summary`)

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

  write_index(args$outdir)
}

main()
