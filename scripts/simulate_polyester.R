#!/usr/bin/env Rscript
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=simulate_polyester
#SBATCH --output=logs/simulate_polyester_%j.log
#SBATCH --error=logs/simulate_polyester_%j.err

## polyester-based data simulation script
suppressPackageStartupMessages({
  # Ensure required packages are available (install missing deps on the fly)
  install_if_missing <- function(pkg, bioc = FALSE) {
    if (requireNamespace(pkg, quietly = TRUE)) return(invisible(TRUE))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }

  install_if_missing("optparse", bioc = FALSE)
  install_if_missing("Biostrings", bioc = TRUE)
  install_if_missing("polyester", bioc = TRUE)
  library(optparse)
  library(Biostrings)
  library(polyester)
})

option_list <- list(
  make_option(c("-t", "--transcripts"), type = "character", help = "Transcript FASTA file (required)"),
  make_option(c("-o", "--outdir"), type = "character", default = "simulated", help = "Output directory"),
  make_option(c("-n", "--n_samples"), type = "integer", default = 1, help = "Number of samples to simulate"),
  make_option(c("-r", "--reads_per_sample"), type = "integer", default = 1e6, help = "Total reads per sample"),
  make_option(c("-l", "--read_length"), type = "integer", default = 100, help = "Read length"),
  make_option(c("--paired"), action = "store_true", default = TRUE, help = "Generate paired-end reads (default: TRUE)"),
  make_option(c("--seed"), type = "integer", default = 1234, help = "Random seed (default: %default)"),
  make_option(c("--log"), type = "character", default = NA, help = "Optional log file path")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$transcripts) || identical(opt$transcripts, "")) {
  stop("--transcripts is required")
}

if (!file.exists(opt$transcripts)) stop(paste0("Transcripts file not found: ", opt$transcripts))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

if (!is.na(opt$log)) {
  con <- file(opt$log, open = "wt")
  sink(con, type = "output")
  sink(con, type = "message")
}

set.seed(opt$seed)

tx <- Biostrings::readDNAStringSet(opt$transcripts)
if (length(tx) == 0) stop("No transcripts found in FASTA file")

# Distribute reads evenly across transcripts as a simple default
k <- length(tx)
reads_per_transcript <- matrix(
  floor(opt$reads_per_sample / k),
  nrow = k, ncol = opt$n_samples
)
# Keep expression constant across samples unless user provides custom design later
fold_changes <- matrix(1, nrow = k, ncol = opt$n_samples)

message("Simulating ", opt$n_samples, " sample(s) with ~", opt$reads_per_sample, " reads each")

simulate_experiment(
  fasta = opt$transcripts,
  readmat = reads_per_transcript,
  fold_changes = fold_changes,
  outdir = opt$outdir,
  readlen = opt$read_length,
  paired = opt$paired,
  gzip = TRUE
)

message("Done! Simulated reads are in ", opt$outdir)

if (!is.na(opt$log)) {
  sink(type = "message")
  sink(type = "output")
  close(con)
}
