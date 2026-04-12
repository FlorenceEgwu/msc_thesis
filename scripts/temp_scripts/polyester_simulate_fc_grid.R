#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
})

# -----------------------------
# Settings
# -----------------------------
READLEN <- 100
COVERAGE_LEVELS <- c(5, 10, 20)     # C(5,10,20) #change in coverage more important.
REPLICATE_LEVELS <- c(5)    # rpt x5 x10 x20 #rpt fixed at 5
SEED_BASE <- 12345

# FC groups (Cond1:Cond2)
FC_GROUPS <- list(
  `1_1` = c(1, 1),
  `1_2` = c(1, 2),
  `2_1` = c(2, 1),
  `1_4` = c(1, 4),
  `4_1` = c(4, 1)
)

# Base 300-transcript template:
# 1:1 = 100 transcripts, others = 50 transcripts each
BASE_TEMPLATE <- c(`1_1` = 100, `1_2` = 50, `2_1` = 50, `1_4` = 50, `4_1` = 50)
 
# -----------------------------
# Helpers
# -----------------------------
make_group_sizes <- function(n_transcripts) {
  if (n_transcripts %% 300 != 0) {
    stop("Transcript count must be a multiple of 300 to match 300-template FC design.")
  }
  mult <- n_transcripts / 300
  BASE_TEMPLATE * mult
}

make_fc_matrix <- function(group_sizes) {
  n <- sum(group_sizes)
  fc <- matrix(NA_real_, nrow = n, ncol = 2)
  colnames(fc) <- c("Cond1", "Cond2")

  idx <- 1
  for (g in names(group_sizes)) {
    k <- as.integer(group_sizes[[g]])
    fc[idx:(idx + k - 1), ] <- matrix(FC_GROUPS[[g]], nrow = k, ncol = 2, byrow = TRUE)
    idx <- idx + k
  }
  fc
}

make_group_labels <- function(group_sizes) {
  rep(names(group_sizes), times = as.integer(group_sizes))
}

# reads_per_transcript = coverage * transcript_length / read_length
reads_per_transcript_from_coverage <- function(tx_len, coverage, readlen) {
  rpt <- floor((coverage * tx_len) / readlen)
  rpt[rpt < 1] <- 1
  rpt
}

write_design_table <- function(outdir, tx_ids, tx_len, groups, fc_mat) {
  design <- data.frame(
    transcript_id = tx_ids,
    length = tx_len,
    group = groups,
    cond1_mult = fc_mat[, 1],
    cond2_mult = fc_mat[, 2],
    stringsAsFactors = FALSE
  )
  write.csv(design, file = file.path(outdir, "transcript_fc_design.csv"), row.names = FALSE)
}

simulate_one <- function(fasta_path, out_root, dataset_name) {
  tx <- readDNAStringSet(fasta_path)
  n_tx <- length(tx)

  group_sizes <- make_group_sizes(n_tx)
  fc_mat <- make_fc_matrix(group_sizes)
  groups <- make_group_labels(group_sizes)

  #To ensure ordering is stable and reproducible (use FASTA order)
  tx_ids <- names(tx)
  tx_len <- width(tx)

  # Safety
  if (length(groups) != n_tx) stop("Internal error: groups != number of transcripts.")

  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

  for (cov in COVERAGE_LEVELS) {
    for (rpt in REPLICATE_LEVELS) {

      set.seed(SEED_BASE + cov * 100 + rpt)

      outdir <- file.path(
        out_root,
        dataset_name,
        paste0("COV_", cov),
        paste0("RPT_", rpt)
      )
      dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

      # Use exact FASTA (copy into run folder for provenance)
      fasta_used <- file.path(outdir, "transcripts_used.fa")
      writeXStringSet(tx, filepath = fasta_used, format = "fasta")

      # reads per transcript formula
      reads_per_tx <- reads_per_transcript_from_coverage(tx_len, cov, READLEN)

      # write design metadata
      write_design_table(outdir, tx_ids, tx_len, groups, fc_mat)

      # simulate
      simulate_experiment(
        fasta = fasta_used,
        reads_per_transcript = reads_per_tx,
        num_reps = c(rpt, rpt),          # two conditions, same reps each
        fold_changes = fc_mat,
        readlen = READLEN,
        paired = TRUE,                   # set FALSE if single-end is needed
        strand_specific = FALSE,         # set TRUE if needed
        outdir = outdir
      )

      message(sprintf(
        "[%s] DONE cov=%d rpt=%d -> %s",
        dataset_name, cov, rpt, outdir
      ))
    }
  }
}

# -----------------------------
# Run both datasets
# -----------------------------
# Paths you uploaded (as mounted in this environment)
FASTA_DATASET1 <- "data/input/fasta_polyester/fasta_prep/transcripts_1tx.fa"  # 900 transcripts
FASTA_DATASET2 <- "data/input/fasta_polyester/fasta_prep/transcripts_2tx.fa"  # 1800 transcripts

OUT_ROOT <- "data/sim/fc_grid_design"

simulate_one(FASTA_DATASET1, OUT_ROOT, "dataset1_900tx")
simulate_one(FASTA_DATASET2, OUT_ROOT, "dataset2_1800tx")

message("ALL SIMULATIONS COMPLETE.")