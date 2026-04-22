#!/usr/bin/env Rscript

# Build two datasets with strict transcript/exon composition and non-overlapping genes:
# dataset 1: 900 genes with exactly 1 transcript (300 <=2 exons, 600 >=3 exons)
# dataset 2: 900 genes with exactly 2 transcripts (300 genes where both tx <=2 exons,
#            600 genes where both tx >=3 exons)

suppressPackageStartupMessages({
  required_packages <- c(
    "rtracklayer", "GenomicFeatures", "GenomicRanges", "Biostrings", "Rsamtools",
    "dplyr", "data.table"
  )
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(lapply(required_packages, library, character.only = TRUE))
})

# ----------------------------
# Configuration
# ----------------------------
gtf_path <- "reference/genome/Mus_musculus.GRCm39.115.gtf"
cdna_fa <- "reference/genome/Mus_musculus.GRCm39.cdna.all.fa.gz"
out_dir <- "data/input/selected_transcripts"
protein_coding_only <- TRUE
seed <- 20260220L

n_total <- 900L
n_le2 <- 300L
n_ge3 <- 600L

if (n_le2 + n_ge3 != n_total) {
  stop("n_le2 + n_ge3 must equal n_total.", call. = FALSE)
}

set.seed(seed)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading annotation...\n")
txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")
genes_gr <- genes(txdb)
exons_by_tx <- exonsBy(txdb, by = "tx")

strip_ensembl_version <- function(x) sub("\\..*$", "", x)

cat("Indexing Ensembl cDNA FASTA...\n")
cdna_index <- fasta.index(cdna_fa, seqtype = "DNA")
if (nrow(cdna_index) == 0L) {
  stop("No records found in Ensembl cDNA FASTA: ", cdna_fa, call. = FALSE)
}
cdna_header_ids <- sub(" .*", "", cdna_index$desc)
cdna_header_ids_base <- strip_ensembl_version(cdna_header_ids)
cdna_row_by_id <- setNames(seq_len(nrow(cdna_index)), cdna_header_ids)
cdna_row_by_id <- cdna_row_by_id[!duplicated(names(cdna_row_by_id))]
cdna_row_by_base_id <- setNames(seq_len(nrow(cdna_index)), cdna_header_ids_base)
cdna_row_by_base_id <- cdna_row_by_base_id[!duplicated(names(cdna_row_by_base_id))]

tx_map <- AnnotationDbi::select(
  txdb,
  keys = keys(txdb, keytype = "GENEID"),
  columns = c("TXID", "TXNAME"),
  keytype = "GENEID"
)

tx_tbl <- data.frame(
  gene_id = as.character(tx_map$GENEID),
  tx_id = as.character(tx_map$TXID),
  tx_name = as.character(tx_map$TXNAME),
  stringsAsFactors = FALSE
) %>%
  distinct(gene_id, tx_id, tx_name)

exon_counts_tbl <- data.frame(
  tx_id = as.character(names(exons_by_tx)),
  n_exons = as.integer(lengths(exons_by_tx)),
  stringsAsFactors = FALSE
)

tx_tbl <- tx_tbl %>% left_join(exon_counts_tbl, by = "tx_id")
if (any(is.na(tx_tbl$n_exons))) {
  stop("Could not determine exon counts for some transcripts.", call. = FALSE)
}

gene_tbl <- data.frame(
  gene_id = names(genes_gr),
  seqnames = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  strand = as.character(strand(genes_gr)),
  stringsAsFactors = FALSE
)

cat("Loading gene metadata...\n")
gene_rows <- tryCatch(
  as.data.frame(import(gtf_path, format = "gtf", feature.type = "gene")),
  error = function(e) {
    raw_tbl <- as.data.frame(import(gtf_path))
    raw_tbl[as.character(raw_tbl$type) == "gene", , drop = FALSE]
  }
)

col_or_na <- function(df, col) {
  if (col %in% names(df)) df[[col]] else rep(NA_character_, nrow(df))
}

gene_meta <- data.frame(
  gene_id = coalesce(
    col_or_na(gene_rows, "gene_id"),
    col_or_na(gene_rows, "ID"),
    col_or_na(gene_rows, "gene")
  ),
  gene_name = coalesce(
    col_or_na(gene_rows, "gene_name"),
    col_or_na(gene_rows, "Name")
  ),
  gene_biotype = coalesce(
    col_or_na(gene_rows, "gene_biotype"),
    col_or_na(gene_rows, "gene_type"),
    col_or_na(gene_rows, "biotype")
  ),
  stringsAsFactors = FALSE
)
gene_meta <- gene_meta[!is.na(gene_meta$gene_id), , drop = FALSE]
gene_meta <- gene_meta[!duplicated(gene_meta$gene_id), , drop = FALSE]

gene_tx_summary <- tx_tbl %>%
  group_by(gene_id) %>%
  summarise(
    n_transcripts = n_distinct(tx_id),
    min_exons_per_tx = min(n_exons),
    max_exons_per_tx = max(n_exons),
    all_exons_le2 = all(n_exons <= 2L),
    all_exons_ge3 = all(n_exons >= 3L),
    .groups = "drop"
  )

combined <- gene_tx_summary %>%
  left_join(gene_tbl, by = "gene_id") %>%
  left_join(gene_meta, by = "gene_id")

if (protein_coding_only && "gene_biotype" %in% colnames(combined)) {
  combined <- combined %>% filter(gene_biotype == "protein_coding")
}

valid_gene_ids <- combined$gene_id
tx_tbl <- tx_tbl %>% filter(gene_id %in% valid_gene_ids)
genes_gr <- genes_gr[intersect(names(genes_gr), valid_gene_ids)]

if (length(genes_gr) == 0L) {
  stop("No genes available after filtering.", call. = FALSE)
}

gene_chr <- setNames(as.character(seqnames(genes_gr)), names(genes_gr))
gene_start <- setNames(as.integer(start(genes_gr)), names(genes_gr))
gene_end <- setNames(as.integer(end(genes_gr)), names(genes_gr))

drop_locked_overlaps <- function(candidate_ids, locked_ids) {
  candidate_ids <- unique(candidate_ids)
  locked_ids <- unique(locked_ids)

  if (length(candidate_ids) == 0L || length(locked_ids) == 0L) {
    return(candidate_ids)
  }

  ov <- findOverlaps(
    genes_gr[candidate_ids],
    genes_gr[locked_ids],
    ignore.strand = TRUE
  )
  if (length(ov) == 0L) {
    return(candidate_ids)
  }

  candidate_ids[-unique(queryHits(ov))]
}

max_non_overlapping_ids <- function(candidate_ids) {
  candidate_ids <- unique(candidate_ids)
  if (length(candidate_ids) == 0L) {
    return(character())
  }

  dt <- data.table::data.table(
    gene_id = candidate_ids,
    chr = gene_chr[candidate_ids],
    start = gene_start[candidate_ids],
    end = gene_end[candidate_ids],
    tie = runif(length(candidate_ids))
  )
  data.table::setorder(dt, chr, end, tie, start)

  selected <- character()
  by_chr <- split(dt, dt$chr, drop = TRUE)
  for (chr_dt in by_chr) {
    last_end <- -Inf
    for (i in seq_len(nrow(chr_dt))) {
      if (chr_dt$start[[i]] > last_end) {
        selected <- c(selected, chr_dt$gene_id[[i]])
        last_end <- chr_dt$end[[i]]
      }
    }
  }

  selected
}

pick_non_overlapping <- function(candidate_ids, n, locked_ids = character(), label = "selection") {
  candidate_ids <- drop_locked_overlaps(candidate_ids, locked_ids)
  if (length(candidate_ids) < n) {
    stop(label, ": need ", n, " candidates, found ", length(candidate_ids), ".", call. = FALSE)
  }

  max_set <- max_non_overlapping_ids(candidate_ids)
  if (length(max_set) < n) {
    stop(label, ": only ", length(max_set), " non-overlapping genes available.", call. = FALSE)
  }

  sample(max_set, n, replace = FALSE)
}

assert_no_overlap <- function(gene_ids, label) {
  gr <- genes_gr[unique(gene_ids)]
  ov <- findOverlaps(gr, gr, ignore.strand = TRUE)
  if (any(queryHits(ov) < subjectHits(ov))) {
    stop(label, " has overlaps.", call. = FALSE)
  }
}

cat("Selecting candidate pools...\n")
candidates_1tx_le2 <- combined %>%
  filter(n_transcripts == 1L, max_exons_per_tx <= 2L) %>%
  pull(gene_id)

candidates_1tx_ge3 <- combined %>%
  filter(n_transcripts == 1L, min_exons_per_tx >= 3L) %>%
  pull(gene_id)

candidates_2tx_le2 <- combined %>%
  filter(n_transcripts == 2L, all_exons_le2) %>%
  pull(gene_id)

candidates_2tx_ge3 <- combined %>%
  filter(n_transcripts == 2L, all_exons_ge3) %>%
  pull(gene_id)

cat("Selecting dataset 1...\n")
selected_1tx_le2 <- pick_non_overlapping(candidates_1tx_le2, n_le2, label = "dataset1_le2")
selected_1tx_ge3 <- pick_non_overlapping(candidates_1tx_ge3, n_ge3, locked_ids = selected_1tx_le2, label = "dataset1_ge3")
selected_1tx <- unique(c(selected_1tx_le2, selected_1tx_ge3))
if (length(selected_1tx) != n_total) stop("dataset1 does not contain exactly 900 genes.", call. = FALSE)
assert_no_overlap(selected_1tx, "dataset1")

cat("Selecting dataset 2...\n")
selected_2tx_le2 <- pick_non_overlapping(candidates_2tx_le2, n_le2, label = "dataset2_le2")
selected_2tx_ge3 <- pick_non_overlapping(candidates_2tx_ge3, n_ge3, locked_ids = selected_2tx_le2, label = "dataset2_ge3")
selected_2tx <- unique(c(selected_2tx_le2, selected_2tx_ge3))
if (length(selected_2tx) != n_total) stop("dataset2 does not contain exactly 900 genes.", call. = FALSE)
assert_no_overlap(selected_2tx, "dataset2")

n_tx_dataset1 <- tx_tbl %>%
  filter(gene_id %in% selected_1tx) %>%
  summarise(n = n_distinct(tx_id)) %>%
  pull(n)

n_tx_dataset2 <- tx_tbl %>%
  filter(gene_id %in% selected_2tx) %>%
  summarise(n = n_distinct(tx_id)) %>%
  pull(n)

if (n_tx_dataset1 != 900L) stop("dataset1 transcript count expected 900, found ", n_tx_dataset1, ".", call. = FALSE)
if (n_tx_dataset2 != 1800L) stop("dataset2 transcript count expected 1800, found ", n_tx_dataset2, ".", call. = FALSE)

extract_tx_fasta_ordered <- function(ordered_tx_names, tx_to_gene_map, out_fa) {
  requested_ids <- ordered_tx_names
  if (length(requested_ids) == 0L) return(invisible(NULL))

  if (anyDuplicated(requested_ids) > 0L) {
    stop("Duplicate transcript IDs found in requested FASTA order.", call. = FALSE)
  }

  exact_hits <- unname(cdna_row_by_id[requested_ids])
  base_hits <- unname(cdna_row_by_base_id[strip_ensembl_version(requested_ids)])
  row_idx <- exact_hits
  row_idx[is.na(row_idx)] <- base_hits[is.na(row_idx)]

  missing_ids <- requested_ids[is.na(row_idx)]
  if (length(missing_ids) > 0L) {
    stop(
      "Could not find ", length(missing_ids), " transcript IDs in Ensembl cDNA FASTA for ",
      out_fa, ". Example IDs: ", paste(utils::head(missing_ids, 5L), collapse = ", "),
      call. = FALSE
    )
  }

  if (any(is.na(tx_to_gene_map[requested_ids]))) {
    stop("Missing transcript->gene mapping for ordered FASTA output.", call. = FALSE)
  }

  seqs <- readDNAStringSet(cdna_index[as.integer(row_idx), , drop = FALSE])
  names(seqs) <- paste0(requested_ids, "|", tx_to_gene_map[requested_ids])
  writeXStringSet(seqs, out_fa)
}

ordered_1tx_genes <- c(sort(selected_1tx_le2), sort(selected_1tx_ge3))
ordered_2tx_genes <- c(sort(selected_2tx_le2), sort(selected_2tx_ge3))

tx_order_1tx <- tx_tbl %>%
  filter(gene_id %in% ordered_1tx_genes) %>%
  distinct(gene_id, tx_id, tx_name) %>%
  mutate(
    tx_id_num = suppressWarnings(as.integer(tx_id)),
    gene_order = match(gene_id, ordered_1tx_genes)
  ) %>%
  arrange(gene_order, tx_id_num, tx_name)

if (nrow(tx_order_1tx) != 900L) {
  stop("Internal ordering error: dataset1 expected 900 transcripts.", call. = FALSE)
}
if (anyDuplicated(tx_order_1tx$gene_id) > 0L) {
  stop("Internal ordering error: dataset1 expected one transcript per gene.", call. = FALSE)
}

tx_order_2tx <- tx_tbl %>%
  filter(gene_id %in% ordered_2tx_genes) %>%
  distinct(gene_id, tx_id, tx_name) %>%
  mutate(
    tx_id_num = suppressWarnings(as.integer(tx_id)),
    gene_order = match(gene_id, ordered_2tx_genes)
  ) %>%
  arrange(gene_order, tx_id_num, tx_name) %>%
  group_by(gene_id) %>%
  mutate(transcript_in_gene = row_number()) %>%
  ungroup()

if (nrow(tx_order_2tx) != 1800L) {
  stop("Internal ordering error: dataset2 expected 1800 transcripts.", call. = FALSE)
}
if (!all(table(tx_order_2tx$gene_id) == 2L)) {
  stop("Internal ordering error: dataset2 expected exactly 2 transcripts per gene.", call. = FALSE)
}

ordered_tx_1tx <- tx_order_1tx$tx_name
ordered_tx_2tx <- c(
  tx_order_2tx %>% filter(transcript_in_gene == 1L) %>% arrange(gene_order) %>% pull(tx_name),
  tx_order_2tx %>% filter(transcript_in_gene == 2L) %>% arrange(gene_order) %>% pull(tx_name)
)

tx_to_gene_1tx <- setNames(tx_order_1tx$gene_id, tx_order_1tx$tx_name)
tx_to_gene_2tx <- setNames(tx_order_2tx$gene_id, tx_order_2tx$tx_name)

dataset_1_tbl <- combined %>%
  filter(gene_id %in% ordered_1tx_genes) %>%
  mutate(dataset = "dataset_1", exon_group = ifelse(gene_id %in% selected_1tx_le2, "le2_exons", "ge3_exons")) %>%
  mutate(gene_order = match(gene_id, ordered_1tx_genes)) %>%
  arrange(gene_order) %>%
  select(-gene_order)

dataset_2_tbl <- combined %>%
  filter(gene_id %in% ordered_2tx_genes) %>%
  mutate(dataset = "dataset_2", exon_group = ifelse(gene_id %in% selected_2tx_le2, "le2_exons", "ge3_exons")) %>%
  mutate(gene_order = match(gene_id, ordered_2tx_genes)) %>%
  arrange(gene_order) %>%
  select(-gene_order)

summary_tbl <- data.frame(
  dataset = c("dataset_1", "dataset_2"),
  n_genes = c(nrow(dataset_1_tbl), nrow(dataset_2_tbl)),
  n_transcripts = c(n_tx_dataset1, n_tx_dataset2),
  genes_le2_exons = c(length(selected_1tx_le2), length(selected_2tx_le2)),
  genes_ge3_exons = c(length(selected_1tx_ge3), length(selected_2tx_ge3)),
  stringsAsFactors = FALSE
)

out_1tx_genes <- file.path(out_dir, "genes_1tx.tsv")
out_2tx_genes <- file.path(out_dir, "genes_2tx.tsv")
out_1tx_fa <- file.path(out_dir, "transcripts_1tx.fa")
out_2tx_fa <- file.path(out_dir, "transcripts_2tx.fa")
out_summary <- file.path(out_dir, "gene_selector_summary.tsv")

data.table::fwrite(dataset_1_tbl, out_1tx_genes, sep = "\t")
data.table::fwrite(dataset_2_tbl, out_2tx_genes, sep = "\t")
data.table::fwrite(summary_tbl, out_summary, sep = "\t")

cat("Writing FASTA files...\n") 
extract_tx_fasta_ordered(ordered_tx_1tx, tx_to_gene_1tx, out_1tx_fa)
extract_tx_fasta_ordered(ordered_tx_2tx, tx_to_gene_2tx, out_2tx_fa)

cat("Done.\n")
cat("Wrote:\n")
cat("  ", out_1tx_genes, "\n", sep = "")
cat("  ", out_2tx_genes, "\n", sep = "")
cat("  ", out_1tx_fa, "\n", sep = "")
cat("  ", out_2tx_fa, "\n", sep = "")
cat("  ", out_summary, "\n", sep = "")
