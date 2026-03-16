#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  ensure_pkg <- function(pkg, bioc = FALSE) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      return(invisible(TRUE))
    }

    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }

  ensure_pkg("rtracklayer", bioc = TRUE)
  ensure_pkg("GenomicFeatures", bioc = TRUE)
  ensure_pkg("dplyr")
  ensure_pkg("readr")
  ensure_pkg("tibble")
  ensure_pkg("stringr")

  library(rtracklayer)
  library(GenomicFeatures)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

# ----------------------------
# User configuration
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

pick_default_annotation <- function() {
  candidates <- Sys.glob("reference/Mus_musculus.GRCm39.*.gtf.gz")
  if (length(candidates) == 0) {
    return("reference/Mus_musculus.GRCm39.115.gtf.gz")
  }

  rel <- suppressWarnings(as.integer(str_match(basename(candidates), "GRCm39\\.(\\d+)")[, 2]))
  if (all(is.na(rel))) {
    return(sort(candidates)[[1]])
  }

  candidates[[which.max(rel)]]
}

if (length(args) < 1) {
  cat("Usage:\n")
  cat("  Rscript geneFinder.R <out_prefix> [annotation.gtf|gff3[.gz]] [protein_coding_only TRUE|FALSE]\n\n")
  cat("Examples:\n")
  cat("  Rscript geneFinder.R mouse_grcm39\n")
  cat("  Rscript geneFinder.R mouse_grcm39 Mus_musculus.GRCm39.115.gtf.gz TRUE\n")
  quit(status = 1)
}

out_prefix <- args[[1]]
annotation_path <- ifelse(length(args) >= 2, args[[2]], pick_default_annotation())
protein_coding_only <- ifelse(length(args) >= 3, as.logical(args[[3]]), TRUE)

# ----------------------------
# Helpers
# ----------------------------
is_gtf <- function(path) {
  str_detect(tolower(path), "\\.gtf(\\.gz)?$")
}

is_gff3 <- function(path) {
  str_detect(tolower(path), "\\.gff3?(\\.gz)?$")
}

stop_if_not_supported <- function(path) {
  if (!is_gtf(path) && !is_gff3(path)) {
    stop("Unsupported annotation format. Provide .gtf/.gtf.gz or .gff/.gff3/.gz")
  }
}

extract_gene_and_tx_tables <- function(gtf_or_gff) {
  get_tx_ids <- function(tx_gr) {
    tx_ids <- names(tx_gr)
    if (is.null(tx_ids) || all(is.na(tx_ids)) || all(tx_ids == "")) {
      tx_mcols <- mcols(tx_gr)
      if ("tx_name" %in% colnames(tx_mcols)) {
        tx_ids <- tx_mcols$tx_name
      } else if ("tx_id" %in% colnames(tx_mcols)) {
        tx_ids <- tx_mcols$tx_id
      } else {
        tx_ids <- seq_along(tx_gr)
      }
    }
    as.character(tx_ids)
  }

  # Import as GRanges
  gr <- import(gtf_or_gff)

  # rtracklayer conventions:
  # - GTF: type column often in 'type' (from GRanges), and gene_id/transcript_id in mcols
  # - GFF3: uses ID/Parent or gene_id/transcript_id depending on source
  # We'll build TxDb to standardize extraction.
  txdb <- if (is_gtf(gtf_or_gff)) {
    makeTxDbFromGFF(gtf_or_gff, format = "gtf")
  } else {
    makeTxDbFromGFF(gtf_or_gff, format = "gff3")
  }

  # Genes table (gene_id and coordinates)
  genes_gr <- genes(txdb)
  genes_tbl <- tibble(
    gene_id = names(genes_gr),
    seqnames = as.character(seqnames(genes_gr)),
    start = start(genes_gr),
    end = end(genes_gr),
    strand = as.character(strand(genes_gr))
  )

  raw_tbl <- as_tibble(as.data.frame(gr))

  needed_cols <- c(
    "gene_id", "ID", "gene", "gene_name", "Name",
    "gene_biotype", "gene_type", "biotype",
    "transcript_id", "Parent", "transcript", "transcript_name",
    "exon_id"
  )
  for (col in needed_cols) {
    if (!col %in% names(raw_tbl)) {
      raw_tbl[[col]] <- NA_character_
    }
  }

  # Transcripts with gene mapping (prefer raw annotation IDs)
  tx_tbl <- raw_tbl %>%
    transmute(
      gene_id = coalesce(.data$gene_id, .data$gene),
      transcript_id = coalesce(.data$transcript_id, .data$Parent, .data$transcript, .data$ID)
    ) %>%
    filter(!is.na(gene_id), !is.na(transcript_id)) %>%
    distinct(gene_id, transcript_id)

  if (nrow(tx_tbl) == 0) {
    tx_by_gene <- transcriptsBy(txdb, by = "gene")
    tx_tbl <- bind_rows(lapply(names(tx_by_gene), function(g) {
      tx_gr <- tx_by_gene[[g]]
      tx_ids <- get_tx_ids(tx_gr)
      tibble(
        gene_id = g,
        transcript_id = as.character(tx_ids)
      )
    }))
  }

  exon_counts_tbl <- raw_tbl %>%
    filter(type == "exon") %>%
    transmute(
      transcript_id = coalesce(.data$transcript_id, .data$Parent, .data$transcript, .data$ID),
      exon_id = coalesce(.data$exon_id, .data$ID)
    ) %>%
    filter(!is.na(transcript_id)) %>%
    group_by(transcript_id) %>%
    summarise(n_exons = n(), .groups = "drop")

  if (nrow(exon_counts_tbl) == 0) {
    exons_by_tx <- exonsBy(txdb, by = "tx")
    tx_ids_exons <- get_tx_ids(exons_by_tx)
    exon_counts_tbl <- tibble(
      transcript_id = tx_ids_exons,
      n_exons = lengths(exons_by_tx)
    )
  }

  # Try to build a gene metadata table from rows of type "gene" (GTF) or features with gene-level IDs
  gene_meta <- raw_tbl %>%
    mutate(type = as.character(type)) %>%
    filter(type == "gene") %>%
    transmute(
      gene_id = coalesce(.data$gene_id, .data$ID, .data$gene),
      gene_name = coalesce(.data$gene_name, .data$Name),
      gene_biotype = coalesce(.data$gene_biotype, .data$gene_type, .data$biotype)
    ) %>%
    filter(!is.na(gene_id)) %>%
    distinct(gene_id, .keep_all = TRUE)

  list(
    genes_tbl = genes_tbl,
    tx_tbl = tx_tbl,
    exon_counts_tbl = exon_counts_tbl,
    gene_meta = gene_meta
  )
}

# ----------------------------
# Main
# ----------------------------
stop_if_not_supported(annotation_path)

cat("Reading annotation: ", annotation_path, "\n", sep = "")
cat("protein_coding_only: ", protein_coding_only, "\n", sep = "")

tables <- extract_gene_and_tx_tables(annotation_path)

genes_tbl <- tables$genes_tbl
tx_tbl    <- tables$tx_tbl
exon_counts_tbl <- tables$exon_counts_tbl
gene_meta <- tables$gene_meta

# Count unique transcripts per gene
counts_tbl <- tx_tbl %>%
  distinct(gene_id, transcript_id) %>%
  group_by(gene_id) %>%
  summarise(n_transcripts = n(), .groups = "drop")

tx_exon_tbl <- tx_tbl %>%
  distinct(gene_id, transcript_id) %>%
  left_join(exon_counts_tbl, by = "transcript_id")

safe_stat <- function(x, fn) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  fn(x, na.rm = TRUE)
}

tx_exon_summary <- tx_exon_tbl %>%
  group_by(gene_id) %>%
  summarise(
    min_exons_per_tx = safe_stat(n_exons, min),
    max_exons_per_tx = safe_stat(n_exons, max),
    mean_exons_per_tx = safe_stat(n_exons, mean),
    median_exons_per_tx = safe_stat(n_exons, median),
    .groups = "drop"
  )

# Combine tables
combined <- counts_tbl %>%
  left_join(genes_tbl, by = "gene_id") %>%
  left_join(tx_exon_summary, by = "gene_id") %>%
  left_join(gene_meta, by = "gene_id")

# Optional filter: protein-coding genes only (depends on annotation content)
if (protein_coding_only) {
  if (!("gene_biotype" %in% colnames(combined))) {
    warning("gene_biotype not found in annotation; cannot filter protein_coding_only. Proceeding without filter.")
  } else {
    combined <- combined %>%
      filter(!is.na(gene_biotype)) %>%
      filter(gene_biotype == "protein_coding")
  }
}

genes_1tx <- combined %>% filter(n_transcripts == 1) %>% arrange(gene_id)
genes_2tx <- combined %>% filter(n_transcripts == 2) %>% arrange(gene_id)

genes_1tx_exon_summary <- tx_exon_tbl %>%
  filter(gene_id %in% genes_1tx$gene_id) %>%
  group_by(n_exons) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    n_transcripts = n(),
    .groups = "drop"
  ) %>%
  arrange(n_exons)

genes_2tx_exon_summary <- tx_exon_tbl %>%
  filter(gene_id %in% genes_2tx$gene_id) %>%
  group_by(n_exons) %>%
  summarise(
    n_genes = n_distinct(gene_id),
    n_transcripts = n(),
    .groups = "drop"
  ) %>%
  arrange(n_exons)

# Write outputs
out_all  <- paste0(out_prefix, "_gene_transcript_counts.tsv")
out_1tx  <- paste0(out_prefix, "_genes_1_transcript.tsv")
out_2tx  <- paste0(out_prefix, "_genes_2_transcripts.tsv")
out_tx_exons <- paste0(out_prefix, "_transcript_exon_counts.tsv")
out_1tx_exon_summary <- paste0(out_prefix, "_genes_1_transcript_exon_summary.tsv")
out_2tx_exon_summary <- paste0(out_prefix, "_genes_2_transcripts_exon_summary.tsv")
out_sum  <- paste0(out_prefix, "_summary.txt")

write_tsv(combined, out_all)
write_tsv(genes_1tx, out_1tx)
write_tsv(genes_2tx, out_2tx)
write_tsv(tx_exon_tbl, out_tx_exons)
write_tsv(genes_1tx_exon_summary, out_1tx_exon_summary)
write_tsv(genes_2tx_exon_summary, out_2tx_exon_summary)

summary_lines <- c(
  paste0("Annotation: ", annotation_path),
  paste0("protein_coding_only: ", protein_coding_only),
  paste0("Total genes considered: ", nrow(combined)),
  paste0("Genes with 1 transcript: ", nrow(genes_1tx)),
  paste0("Genes with 2 transcripts: ", nrow(genes_2tx))
)
writeLines(summary_lines, out_sum)

cat("Done.\n")
cat("Wrote:\n")
cat("  ", out_all, "\n", sep = "")
cat("  ", out_1tx, "\n", sep = "")
cat("  ", out_2tx, "\n", sep = "")
cat("  ", out_tx_exons, "\n", sep = "")
cat("  ", out_1tx_exon_summary, "\n", sep = "")
cat("  ", out_2tx_exon_summary, "\n", sep = "")
cat("  ", out_sum, "\n", sep = "")
