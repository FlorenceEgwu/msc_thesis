#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(rtracklayer)
  library(dplyr)
  library(Biostrings)
  library(stringr)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--gtf", type = "character"),
  make_option("--complex-exon-threshold", type = "integer", default = 2),
  make_option("--out", type = "character")
)))

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}
read_selected_transcript_ids <- function() {
  fasta_paths <- c(
    "data/input/selected_transcripts/transcripts_1tx.fa",
    "data/input/selected_transcripts/transcripts_2tx.fa"
  )

  headers <- unlist(
    lapply(fasta_paths, function(path) names(readBStringSet(path, format = "fasta"))),
    use.names = FALSE
  )

  transcript_ids <- str_match(headers, "^([^|[:space:]]+)")[, 2]
  unique(na.omit(transcript_ids))
}

selected_transcript_ids <- read_selected_transcript_ids()

exon_gr <- import(opt$gtf)
exon_gr <- exon_gr[exon_gr$type == "exon"]

exon_tbl <- tibble(
  seqnames = as.character(seqnames(exon_gr)),
  start = start(exon_gr),
  end = end(exon_gr),
  strand = as.character(strand(exon_gr)),
  transcript_id = mcols(exon_gr)$transcript_id,
  gene_id = mcols(exon_gr)$gene_id,
  exon_number = suppressWarnings(as.integer(mcols(exon_gr)$exon_number))
) %>%
  filter(!is.na(transcript_id)) %>%
  {
    if (length(selected_transcript_ids) > 0) {
      filter(., transcript_id %in% selected_transcript_ids)
    } else {
      .
    }
  } %>%
  arrange(transcript_id, exon_number, start)

tx_tbl <- exon_tbl %>%
  count(transcript_id, gene_id, name = "exon_count") %>%
  mutate(transcript_type = if_else(exon_count <= opt$`complex-exon-threshold`, "easy", "complex"))

exon_tbl %>%
  left_join(tx_tbl, by = c("transcript_id", "gene_id")) %>%
  write_tsv(opt$out)
