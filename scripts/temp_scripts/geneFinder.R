#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

if (length(args) < 2) {
  cat("Usage:\n")
  cat("  Rscript extract_1x2x_transcript_genes.R <annotation.gtf|gff3[.gz]> <out_prefix> [protein_coding_only TRUE|FALSE]\n\n")
  cat("Example:\n")
  cat("  Rscript extract_1x2x_transcript_genes.R Mus_musculus.GRCm39.113.gtf.gz mouse_grcm39 TRUE\n")
  quit(status = 1)
}

annotation_path <- args[[1]]
out_prefix      <- args[[2]]
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
  # Import as GRanges
  gr <- import(gtf_or_gff)

  mcols_df <- as_tibble(as.data.frame(mcols(gr)))

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

  # Transcripts with gene mapping
  tx_by_gene <- transcriptsBy(txdb, by = "gene")
  tx_tbl <- bind_rows(lapply(names(tx_by_gene), function(g) {
    tx_gr <- tx_by_gene[[g]]
    tibble(
      gene_id = g,
      transcript_id = names(tx_gr)
    )
  }))

  # Pull optional gene_name / gene_biotype if present in raw import columns
  # Ensembl GTF usually has: gene_name, gene_biotype (or gene_type), transcript_biotype (or transcript_type)
  raw <- as.data.frame(gr)
  raw_tbl <- as_tibble(raw)

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
gene_meta <- tables$gene_meta

# Count unique transcripts per gene
counts_tbl <- tx_tbl %>%
  distinct(gene_id, transcript_id) %>%
  group_by(gene_id) %>%
  summarise(n_transcripts = n(), .groups = "drop")

# Combine tables
combined <- counts_tbl %>%
  left_join(genes_tbl, by = "gene_id") %>%
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

# Write outputs
out_all  <- paste0(out_prefix, "_gene_transcript_counts.tsv")
out_1tx  <- paste0(out_prefix, "_genes_1_transcript.tsv")
out_2tx  <- paste0(out_prefix, "_genes_2_transcripts.tsv")
out_sum  <- paste0(out_prefix, "_summary.txt")

write_tsv(combined, out_all)
write_tsv(genes_1tx, out_1tx)
write_tsv(genes_2tx, out_2tx)

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
cat("  ", out_sum, "\n", sep = "")
