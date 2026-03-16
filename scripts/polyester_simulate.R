#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
})

# ============================================================
# Polyester simulation script aligned to resources/settings4sim1.html
# ------------------------------------------------------------
# Key fixes relative to the previous version:
# 1) num_reps is fixed at c(5, 5), as described in the HTML.
# 2) The tables in the HTML are used to ASSIGN transcripts to
#    (complexity x fold-change x coverage) bins.
#    They are NOT extra simulation loops over num_reps.
# 3) reads_per_transcript uses round(cov * width(fasta) / read_length).
# 4) Dataset 2 is handled asymmetrically:
#    - transcripts 1:900 (first transcript in each gene) are all 1:1
#    - transcripts 901:1800 (second transcript in each gene) follow Dataset 1 rules
# 5) Design files are written as TSV:
#    - transcript_design.tsv
#    - gene_design.tsv
#    - sample_design.tsv
#    - count_matrix.tsv
# ============================================================

# -----------------------------
# Settings
# -----------------------------
READ_LENGTH <- 100L
NUM_REPS <- c(5L, 5L) # Cond1 and Cond2 replicates fixed at 5 each
PAIRED <- TRUE
STRAND_SPECIFIC <- FALSE
SEED_BASE <- 12345L

FASTA_DATASET1 <- "data/input/fasta_polyester/fasta_prep_cdna/transcripts_1tx.fa"  # 900 transcripts
FASTA_DATASET2 <- "data/input/fasta_polyester/fasta_prep_cdna/transcripts_2tx.fa"  # 1800 transcripts
OUT_ROOT <- "data/sim/settings4sim1_design"

# FC groups (Cond1:Cond2)
FC_GROUPS <- list(
  "1_1" = c(1, 1),
  "1_2" = c(1, 2),
  "2_1" = c(2, 1),
  "1_4" = c(1, 4),
  "4_1" = c(4, 1)
)

# -----------------------------
# Helpers
# -----------------------------
stopf <- function(...) stop(sprintf(...), call. = FALSE)

reads_per_transcript_from_coverage <- function(tx_len, coverage, read_length) {
  rpt <- round(coverage * tx_len / read_length)
  rpt[rpt < 1] <- 1L
  as.integer(rpt)
}

fc_label_to_pair <- function(fc_label) {
  FC_GROUPS[[fc_label]]
}

build_sample_design <- function(num_reps) {
  cond_names <- c("Cond1", "Cond2")
  sample_ids <- c(
    sprintf("Cond1_rep%02d", seq_len(num_reps[1])),
    sprintf("Cond2_rep%02d", seq_len(num_reps[2]))
  )
  condition <- c(rep("Cond1", num_reps[1]), rep("Cond2", num_reps[2]))
  replicate <- c(seq_len(num_reps[1]), seq_len(num_reps[2]))

  data.frame(
    sample_id = sample_ids,
    condition = factor(condition, levels = cond_names),
    replicate = replicate,
    mate1_fastq = paste0(sample_ids, "_1.fasta"),
    mate2_fastq = if (PAIRED) paste0(sample_ids, "_2.fasta") else NA_character_,
    stringsAsFactors = FALSE
  )
}

make_count_matrix <- function(reads_per_tx, fc_mat, num_reps) {
  sample_design <- build_sample_design(num_reps)
  sample_ids <- sample_design$sample_id

  count_mat <- matrix(
    0L,
    nrow = length(reads_per_tx),
    ncol = length(sample_ids),
    dimnames = list(NULL, sample_ids)
  )

  # Polyester interprets fold_changes as multiplicative factors by condition.
  # This table stores the expected counts passed to each sample in each condition.
  if (num_reps[1] > 0) {
    count_mat[, seq_len(num_reps[1])] <- reads_per_tx * fc_mat[, 1]
  }
  if (num_reps[2] > 0) {
    idx2 <- num_reps[1] + seq_len(num_reps[2])
    count_mat[, idx2] <- reads_per_tx * fc_mat[, 2]
  }

  as.data.frame(count_mat, check.names = FALSE, stringsAsFactors = FALSE)
}

# Infer gene IDs only from ordering described in the settings html, 
# so the we can use the same FASTA files without relying on specific header parsing. 
make_gene_design_dataset1 <- function(tx_ids) {
  n_tx <- length(tx_ids)
  if (n_tx != 900L) stopf("Dataset 1 expects 900 transcripts, got %d.", n_tx)

  data.frame(
    gene_id = sprintf("gene_%03d", seq_len(900L)),
    transcript_id = tx_ids,
    transcript_in_gene = 1L,
    transcript_type = c(rep("easy", 300L), rep("complex", 600L)),
    stringsAsFactors = FALSE
  )
}

make_gene_design_dataset2 <- function(tx_ids) {
  n_tx <- length(tx_ids)
  if (n_tx != 1800L) stopf("Dataset 2 expects 1800 transcripts, got %d.", n_tx)

  gene_id <- sprintf("gene_%03d", rep(seq_len(900L), times = 2L))
  transcript_in_gene <- c(rep(1L, 900L), rep(2L, 900L))
  transcript_type <- c(
    rep("easy", 300L), rep("complex", 600L),
    rep("easy", 300L), rep("complex", 600L)
  )

  data.frame(
    gene_id = gene_id,
    transcript_id = tx_ids,
    transcript_in_gene = transcript_in_gene,
    transcript_type = transcript_type,
    stringsAsFactors = FALSE
  )
}

# Expand a block specification into transcript-level design rows.
# spec columns:
#   transcript_type, fc_label, coverage, n
expand_block_spec <- function(spec_df, expected_n, block_name) {
  rows <- vector("list", nrow(spec_df))
  idx <- 1L

  for (i in seq_len(nrow(spec_df))) {
    row <- spec_df[i, , drop = FALSE]
    n_i <- as.integer(row$n)
    fc_pair <- fc_label_to_pair(row$fc_label)

    rows[[i]] <- data.frame(
      transcript_type = rep(row$transcript_type, n_i),
      fc_label = rep(row$fc_label, n_i),
      cond1_mult = rep(fc_pair[1], n_i),
      cond2_mult = rep(fc_pair[2], n_i),
      coverage = rep(as.integer(row$coverage), n_i),
      stringsAsFactors = FALSE
    )
    idx <- idx + n_i
  }

  out <- do.call(rbind, rows)
  if (nrow(out) != expected_n) {
    stopf(
      "Internal design error for %s: expected %d rows, got %d.",
      block_name, expected_n, nrow(out)
    )
  }
  out
}

# Dataset 1 design (900 transcripts)
# Order required by the HTML:
#   easy x 300, then complex x 600
build_dataset1_design <- function(tx_ids, tx_len) {
  if (length(tx_ids) != 900L) stopf("Dataset 1 expects 900 transcripts.")

  spec_easy <- data.frame(
    transcript_type = "easy",
    fc_label = c("1_1", "1_2", "2_1", "1_4", "4_1"),
    coverage = c(5L, 5L, 5L, 5L, 5L),
    n = c(40L, 15L, 15L, 15L, 15L),
    stringsAsFactors = FALSE
  )
  spec_easy <- do.call(rbind, lapply(c(5L, 10L, 20L), function(cov) {
    transform(spec_easy, coverage = cov)
  }))

  # Reorder by coverage blocks as shown in the table:
  # for each coverage, allocate the five FC groups
  spec_easy <- spec_easy[order(spec_easy$coverage, match(spec_easy$fc_label, c("1_1", "1_2", "2_1", "1_4", "4_1"))), ]

  spec_complex <- data.frame(
    transcript_type = "complex",
    fc_label = c("1_1", "1_2", "2_1", "1_4", "4_1"),
    coverage = c(5L, 5L, 5L, 5L, 5L),
    n = c(80L, 30L, 30L, 30L, 30L),
    stringsAsFactors = FALSE
  )
  spec_complex <- do.call(rbind, lapply(c(5L, 10L, 20L), function(cov) {
    transform(spec_complex, coverage = cov)
  }))
  spec_complex <- spec_complex[order(spec_complex$coverage, match(spec_complex$fc_label, c("1_1", "1_2", "2_1", "1_4", "4_1"))), ]

  design_easy <- expand_block_spec(spec_easy, expected_n = 300L, block_name = "dataset1_easy")
  design_complex <- expand_block_spec(spec_complex, expected_n = 600L, block_name = "dataset1_complex")

  design <- rbind(design_easy, design_complex)
  design$transcript_id <- tx_ids
  design$length <- tx_len
  design$transcript_index <- seq_len(nrow(design))
  design$dataset <- "dataset1"

  design <- design[, c(
    "dataset", "transcript_index", "transcript_id", "length",
    "transcript_type", "coverage", "fc_label", "cond1_mult", "cond2_mult"
  )]

  design
}

# Dataset 2 design (1800 transcripts)
# Order required:
#   1:900  = first transcript in gene: easy x300, complex x600, all FC 1:1
#   901:1800 = second transcript in gene: same rule as Dataset 1
build_dataset2_design <- function(tx_ids, tx_len) {
  if (length(tx_ids) != 1800L) stopf("Dataset 2 expects 1800 transcripts.")

  # First transcript in pair (all FC 1:1)
  spec_first_easy <- data.frame(
    transcript_type = "easy",
    fc_label = "1_1",
    coverage = c(5L, 10L, 20L),
    n = c(100L, 100L, 100L),
    stringsAsFactors = FALSE
  )
  spec_first_complex <- data.frame(
    transcript_type = "complex",
    fc_label = "1_1",
    coverage = c(5L, 10L, 20L),
    n = c(200L, 200L, 200L),
    stringsAsFactors = FALSE
  )

  first_block <- rbind(
    expand_block_spec(spec_first_easy, expected_n = 300L, block_name = "dataset2_first_easy"),
    expand_block_spec(spec_first_complex, expected_n = 600L, block_name = "dataset2_first_complex")
  )

  # Second transcript in pair = Dataset 1 layout
  dummy_ids <- paste0("dummy_", seq_len(900L))
  dummy_len <- rep(100L, 900L)
  second_block <- build_dataset1_design(dummy_ids, dummy_len)
  second_block <- second_block[, c("transcript_type", "coverage", "fc_label", "cond1_mult", "cond2_mult")]

  design <- rbind(first_block, second_block)
  design$transcript_id <- tx_ids
  design$length <- tx_len
  design$transcript_index <- seq_len(nrow(design))
  design$dataset <- "dataset2"
  design$transcript_in_gene <- c(rep(1L, 900L), rep(2L, 900L))

  design <- design[, c(
    "dataset", "transcript_index", "transcript_id", "length",
    "transcript_in_gene", "transcript_type", "coverage",
    "fc_label", "cond1_mult", "cond2_mult"
  )]

  design
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

simulate_dataset <- function(fasta_path, outdir, dataset_name) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tx <- readDNAStringSet(fasta_path)
  tx_ids <- names(tx)
  tx_len <- width(tx)

  if (dataset_name == "dataset1_900tx") {
    design <- build_dataset1_design(tx_ids, tx_len)
    gene_design <- make_gene_design_dataset1(tx_ids)
  } else if (dataset_name == "dataset2_1800tx") {
    design <- build_dataset2_design(tx_ids, tx_len)
    gene_design <- make_gene_design_dataset2(tx_ids)
  } else {
    stopf("Unknown dataset_name: %s", dataset_name)
  }

  if (nrow(design) != length(tx)) {
    stopf(
      "%s design has %d rows but FASTA has %d transcripts.",
      dataset_name, nrow(design), length(tx)
    )
  }

  # Provenance copy of fasta
  fasta_used <- file.path(outdir, "transcripts_used.fa")
  writeXStringSet(tx, filepath = fasta_used, format = "fasta")

  reads_per_tx <- reads_per_transcript_from_coverage(
    tx_len = design$length,
    coverage = design$coverage,
    read_length = READ_LENGTH
  )

  fc_mat <- as.matrix(design[, c("cond1_mult", "cond2_mult")])
  mode(fc_mat) <- "numeric"

  transcript_design <- design
  transcript_design$reads_per_transcript <- reads_per_tx

  sample_design <- build_sample_design(NUM_REPS)
  count_matrix <- make_count_matrix(reads_per_tx, fc_mat, NUM_REPS)
  count_matrix <- cbind(transcript_id = design$transcript_id, count_matrix, stringsAsFactors = FALSE)

  write_tsv(transcript_design, file.path(outdir, "transcript_design.tsv"))
  write_tsv(gene_design, file.path(outdir, "gene_design.tsv"))
  write_tsv(sample_design, file.path(outdir, "sample_design.tsv"))
  write_tsv(count_matrix, file.path(outdir, "count_matrix.tsv"))

  set.seed(SEED_BASE + ifelse(dataset_name == "dataset1_900tx", 1L, 2L))

  simulate_experiment(
    fasta = fasta_used,
    reads_per_transcript = reads_per_tx,
    num_reps = NUM_REPS,
    fold_changes = fc_mat,
    readlen = READ_LENGTH,
    paired = PAIRED,
    strand_specific = STRAND_SPECIFIC,
    outdir = outdir
  )

  message(sprintf("[%s] complete -> %s", dataset_name, outdir))
}

# -----------------------------
# Run
# -----------------------------
simulate_dataset(
  fasta_path = FASTA_DATASET1,
  outdir = file.path(OUT_ROOT, "dataset1_900tx"),
  dataset_name = "dataset1_900tx"
)

simulate_dataset(
  fasta_path = FASTA_DATASET2,
  outdir = file.path(OUT_ROOT, "dataset2_1800tx"),
  dataset_name = "dataset2_1800tx"
)

message("ALL SIMULATIONS COMPLETE.")
