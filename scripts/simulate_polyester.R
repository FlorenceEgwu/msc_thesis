#!/usr/bin/env Rscript
## Robust polyester-based simulation script
suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) stop("Please install R package 'optparse'")
  if (!requireNamespace("Biostrings", quietly = TRUE)) stop("Please install Biostrings")
  if (!requireNamespace("polyester", quietly = TRUE)) stop("Please install polyester")
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

message("Simulating ", opt$n_samples, " sample(s) with ~", opt$reads_per_sample, " reads each")

simulate_experiment(
  fasta = opt$transcripts,
  readmat = reads_per_transcript,
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
