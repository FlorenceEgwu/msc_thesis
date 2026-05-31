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
  make_option("--as-event-tables", type = "character", default = ""),
  make_option("--out", type = "character")
)))

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}
read_selected_transcript_ids <- function() {
  fasta_paths <- c(
    "data/input/selected_transcripts/transcripts_1tx.fa",
    "data/input/selected_transcripts/transcripts_2tx.fa",
    "data/input/selected_transcripts/transcripts_3tx.fa"
  )
  fasta_paths <- fasta_paths[file.exists(fasta_paths)]
  if (length(fasta_paths) == 0) {
    return(character())
  }

  headers <- unlist(
    lapply(fasta_paths, function(path) names(readBStringSet(path, format = "fasta"))),
    use.names = FALSE
  )

  transcript_ids <- str_match(headers, "^([^|[:space:]]+)")[, 2]
  unique(na.omit(transcript_ids))
}

selected_transcript_ids <- read_selected_transcript_ids()

read_as_event_tables <- function(paths_arg) {
  paths_arg <- ifelse(is.null(paths_arg) || is.na(paths_arg), "", paths_arg)
  paths <- strsplit(paths_arg, ",", fixed = TRUE)[[1]]
  paths <- paths[nzchar(paths) & file.exists(paths)]
  if (length(paths) == 0) {
    return(tibble(gene_id = character(), transcript_id = character(),
                  suppa2_as_event_type = character(), suppa2_event_ids = character(),
                  as_event_source = character()))
  }

  bind_rows(lapply(paths, read.delim, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)) %>%
    filter(!is.na(transcript_id), transcript_id != "") %>%
    mutate(
      gene_id = as.character(gene_id),
      transcript_id = as.character(transcript_id),
      as_event_type = coalesce(na_if(as_event_type, ""), "unknown"),
      event_ids = coalesce(event_ids, ""),
      as_event_source = coalesce(na_if(as_event_source, ""), "SUPPA2")
    ) %>%
    group_by(gene_id, transcript_id) %>%
    summarise(
      suppa2_as_event_type = paste(sort(unique(as_event_type)), collapse = ";"),
      suppa2_event_ids = paste(sort(unique(event_ids[event_ids != ""])), collapse = ";"),
      as_event_source = paste(sort(unique(as_event_source)), collapse = ";"),
      .groups = "drop"
    )
}

intron_key <- function(starts, ends, strand) {
  if (length(starts) < 2) {
    return("")
  }
  starts <- as.integer(starts)
  ends <- as.integer(ends)
  if (strand == "-") {
    starts <- rev(starts)
    ends <- rev(ends)
  }
  paste0(head(ends, -1) + 1L, "-", tail(starts, -1) - 1L, collapse = ";")
}

classify_pair_event <- function(tx_a, tx_b) {
  if (nrow(tx_a) == 0 || nrow(tx_b) == 0) {
    return("unknown")
  }
  strand <- tx_a$strand[1]
  introns_a <- intron_key(tx_a$start, tx_a$end, strand)
  introns_b <- intron_key(tx_b$start, tx_b$end, strand)
  if (identical(introns_a, introns_b)) {
    return("same_intron_chain")
  }

  spans_other_intron <- function(exons, other) {
    if (nrow(other) < 2) {
      return(FALSE)
    }
    intron_starts <- head(other$end, -1) + 1L
    intron_ends <- tail(other$start, -1) - 1L
    any(vapply(seq_along(intron_starts), function(i) {
      any(exons$start <= intron_starts[i] & exons$end >= intron_ends[i])
    }, logical(1)))
  }
  if (spans_other_intron(tx_a, tx_b) || spans_other_intron(tx_b, tx_a)) {
    return("RI")
  }

  same_terminal_bounds <- min(tx_a$start) == min(tx_b$start) &&
    max(tx_a$end) == max(tx_b$end)
  if (same_terminal_bounds && nrow(tx_a) != nrow(tx_b)) {
    return("SE")
  }

  if (nrow(tx_a) == nrow(tx_b)) {
    return("A5SS_A3SS_or_MXE")
  }
  "complex_or_other"
}

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

heuristic_as_tbl <- exon_tbl %>%
  group_by(gene_id) %>%
  group_modify(function(.x, .y) {
    tx_ids <- unique(.x$transcript_id)
    if (length(tx_ids) < 2) {
      return(tibble(transcript_id = tx_ids, heuristic_as_event_type = "single_transcript"))
    }
    if (length(tx_ids) > 2) {
      return(tibble(transcript_id = tx_ids, heuristic_as_event_type = "multi_transcript_gene"))
    }
    tx_a <- .x %>% filter(transcript_id == tx_ids[1]) %>% arrange(start, end)
    tx_b <- .x %>% filter(transcript_id == tx_ids[2]) %>% arrange(start, end)
    event_type <- classify_pair_event(tx_a, tx_b)
    tibble(transcript_id = tx_ids, heuristic_as_event_type = event_type)
  }) %>%
  ungroup()

external_as_tbl <- read_as_event_tables(opt$`as-event-tables`)

exon_tbl %>%
  left_join(tx_tbl, by = c("transcript_id", "gene_id")) %>%
  left_join(heuristic_as_tbl, by = c("gene_id", "transcript_id")) %>%
  left_join(external_as_tbl, by = c("gene_id", "transcript_id")) %>%
  mutate(
    as_event_type = if_else(
      !is.na(suppa2_as_event_type) & suppa2_as_event_type != "",
      suppa2_as_event_type,
      heuristic_as_event_type
    ),
    as_event_source = if_else(
      !is.na(suppa2_as_event_type) & suppa2_as_event_type != "",
      as_event_source,
      "heuristic"
    ),
    suppa2_event_ids = if_else(is.na(suppa2_event_ids), "", suppa2_event_ids)
  ) %>%
  write_tsv(opt$out)
