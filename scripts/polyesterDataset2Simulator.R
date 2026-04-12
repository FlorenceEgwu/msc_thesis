#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
  library(dplyr)
})

# Settings and  Configuration
read_length <- 100
num_reps <- c(5, 5) # number of replicates per condition
coverage <- c(5, 10, 20) # target coverage levels for transcripts
paired <- TRUE
set.seed(12345)

fasta_file <- "data/input/selected_transcripts/transcripts_2tx.fa"
tx <- readDNAStringSet(fasta_file)
dataset_name <- "dataset2"
outdir <- "data/input/sim/polyester_design/dataset2"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# helper function for writing files
write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "	", row.names = FALSE, quote = FALSE)
}

# prepare FASTA header parsing to extract transcript and gene IDs
fasta_header_parts <- strsplit(sub(" .*", "", names(tx)), "\\|")
transcript_id <- vapply(fasta_header_parts, `[`, character(1), 1L)
gene_id <- vapply(fasta_header_parts, `[`, character(1), 2L)
first_gene_ids <- gene_id[seq_len(900)]
second_gene_ids <- gene_id[900 + seq_len(900)]

# validate for length & uniqueness of transcripts
if (length(transcript_id) != 1800) {
  stop(sprintf("Dataset 2 expects 1800 transcripts."))
}
if (anyDuplicated(transcript_id) > 0) {
  dup_id <- transcript_id[duplicated(transcript_id)][1]
  stop(sprintf("Duplicate transcript_id found in FASTA: %s", dup_id))
}
# validate genes
if (anyDuplicated(first_gene_ids) > 0 || anyDuplicated(second_gene_ids) > 0) {
  stop("FASTA must contain unique gene IDs within each 900-transcript block.")
}
if (!identical(first_gene_ids, second_gene_ids)) {
  stop("FASTA gene order must match between both 900-transcript blocks.")
}

# FC groups (Cond1:Cond2)
fc_groups <- list(
  "1_1" = c(1, 1),
  "1_2" = c(1, 2),
  "2_1" = c(2, 1),
  "1_4" = c(1, 4),
  "4_1" = c(4, 1)
)
get_fc_pair <- function(fc_label) {
  fc_groups[[fc_label]]
}
# Helper for expanding each tx block spec into tx-level design rows.
expand_block <- function(block_df, expected_n, block_name) {
  out <- block_df[
    rep(seq_len(nrow(block_df)), block_df$n),
    c("transcript_type", "fc_label", "coverage")
  ]
  fc_pairs <- t(vapply(out$fc_label, get_fc_pair, numeric(2)))
  out$cond1_mult <- fc_pairs[, 1]
  out$cond2_mult <- fc_pairs[, 2]
  out$coverage <- as.integer(out$coverage)
  rownames(out) <- NULL # reset row names after expansion
  if (nrow(out) != expected_n) {
    stop(sprintf(
      "Internal design error for %s: expected %d rows, got %d.",
      block_name, expected_n, nrow(out)
    ))
  }
  out
}

# First transcript group in gene (all FC 1:1)
group1_easy <- data.frame(
  transcript_type = "easy",
  fc_label = "1_1",
  coverage = coverage,
  n = c(100, 100, 100)
)
block1_easy <- expand_block(group1_easy, 300, "block1_easy")

group1_complex <- data.frame(
  transcript_type = "complex",
  fc_label = "1_1",
  coverage = coverage,
  n = c(200, 200, 200)
)
block1_complex <- expand_block(group1_complex, 600, "block1_complex")

# First transcript group in gene (multiple FCs)
fc_label <- c("1_1", "1_2", "2_1", "1_4", "4_1")
group2_easy <- data.frame(
  transcript_type = "easy",
  fc_label = rep(fc_label, 3),
  coverage = rep(coverage, 5),
  n = rep(c(40L, 15L, 15L, 15L, 15L), 3)
) %>% arrange(fc_label, coverage)
block2_easy <- expand_block(group2_easy, 300, "block2_easy")

group2_complex <- data.frame(
  transcript_type = "complex",
  fc_label = rep(fc_label, 3),
  coverage = rep(coverage, 5),
  n = rep(c(80L, 30L, 30L, 30L, 30L), 3)
) %>% arrange(fc_label, coverage)
block2_complex <- expand_block(group2_complex, 600, "block2_complex")

# transcript design
tx_design <- rbind(block1_easy, block1_complex, block2_easy, block2_complex)
tx_design$transcript_id <- transcript_id
tx_design$length <- width(tx)
tx_design$transcript_index <- seq_len(nrow(tx_design))
tx_design$dataset <- dataset_name
tx_design$transcript_gene_id <- gene_id
tx_design$fasta_header <- names(tx)
# calc rpt based on coverage and transcript length
reads_per_tx <- round(tx_design$coverage * tx_design$length / read_length)
tx_design$reads_per_tx <- reads_per_tx

# confirm transcript design completeness, order and write to file
if (nrow(tx_design) != length(tx)) {
  stop("%s design has %d rows but FASTA has %d transcripts.",
       dataset_name, nrow(tx_design), length(tx))
}
tx_design <- tx_design[, c(
  "dataset", "transcript_index", "transcript_id", "length", "reads_per_tx",
  "transcript_type", "coverage", "fc_label", "cond1_mult", "cond2_mult",
  "transcript_gene_id", "fasta_header"
)]
write_tsv(tx_design, file.path(outdir, "transcript_design.tsv"))

# Expected sample design
sample_ids <- c(
  sprintf("Cond1_rep%02d", seq_len(num_reps[1])),
  sprintf("Cond2_rep%02d", seq_len(num_reps[2]))
)
sample_design <- data.frame(
  sample_id = sample_ids,
  condition = factor(
    c(rep("Cond1", num_reps[1]), rep("Cond2", num_reps[2])),
    levels = c("Cond1", "Cond2")
  ),
  replicate = c(seq_len(num_reps[1]), seq_len(num_reps[2])),
  mate1_fasta = paste0(sample_ids, "_1.fasta"),
  mate2_fasta = if (paired) paste0(sample_ids, "_2.fasta") else NA_character_,
  stringsAsFactors = FALSE
)
write_tsv(sample_design, file.path(outdir, "sample_design.tsv"))

# expected count matrix
fc_mat <- as.matrix(tx_design[, c("cond1_mult", "cond2_mult")])
mode(fc_mat) <- "numeric"
count_matrix <- matrix(0L,
  nrow = length(tx_design$reads_per_tx),
  ncol = length(sample_ids),
  dimnames = list(NULL, sample_ids)
)
count_matrix[, seq_len(num_reps[1])] <- tx_design$reads_per_tx * fc_mat[, 1]
idx2 <- num_reps[1] + seq_len(num_reps[2])
count_matrix[, idx2] <- tx_design$reads_per_tx * fc_mat[, 2]
count_matrix <- as.data.frame(count_matrix, check.names = FALSE)
count_matrix <- cbind(transcript_id = tx_design$transcript_id, count_matrix)
write_tsv(count_matrix, file.path(outdir, "count_matrix.tsv"))

# Run simulation
message(sprintf("starting simulation for %s", dataset_name))
simulate_experiment(fasta = fasta_file,
  reads_per_transcript = tx_design$reads_per_tx,
  num_reps = num_reps, fold_changes = fc_mat,
  readlen = read_length, paired = paired,
  strand_specific = FALSE, outdir = outdir
)
message(sprintf("dataset simulation complete, check -> %s", outdir))
