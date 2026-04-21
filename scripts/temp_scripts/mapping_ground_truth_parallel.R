#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicAlignments)
  library(Rsamtools)
  library(data.table)
  library(stringr)
  library(Biostrings)
  library(parallel)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--bam", type = "character"),
  make_option("--gtf-table", type = "character"),
  make_option("--truth-read1", type = "character"),
  make_option("--truth-read2", type = "character"),
  make_option("--sample", type = "character"),
  make_option("--mapper", type = "character"),
  make_option("--dataset", type = "character"),
  make_option("--run", type = "character"),
  make_option("--threads", type = "integer", default = 1),
  make_option("--out-coordinates", type = "character"),
  make_option("--out-standard", type = "character"),
  make_option("--out-ground", type = "character"),
  make_option("--out-stratified", type = "character")
)))

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), file = path, sep = "\t", quote = FALSE, na = "NA")
}

pmap_rows <- function(n, fun, threads) {
  if (n == 0L) {
    return(list())
  }
  if (.Platform$OS.type == "windows" || threads <= 1L) {
    return(lapply(seq_len(n), fun))
  }
  mclapply(seq_len(n), fun, mc.cores = threads)
}

parse_truth_headers <- function(path) {
  reads <- readBStringSet(path, format = "fasta")
  headers <- names(reads)
  widths <- width(reads)

  mate1_range <- str_match(headers, "mate1:([0-9]+)-([0-9]+)")
  mate2_range <- str_match(headers, "mate2:([0-9]+)-([0-9]+)")
  mate1_start <- as.integer(mate1_range[, 2])
  mate1_end <- as.integer(mate1_range[, 3])
  mate2_start <- as.integer(mate2_range[, 2])
  mate2_end <- as.integer(mate2_range[, 3])

  mate1_start_only <- as.integer(str_match(headers, "mate1Start:([0-9]+)")[, 2])
  mate2_start_only <- as.integer(str_match(headers, "mate2Start:([0-9]+)")[, 2])

  fill1 <- is.na(mate1_start) & !is.na(mate1_start_only)
  fill2 <- is.na(mate2_start) & !is.na(mate2_start_only)
  mate1_start[fill1] <- mate1_start_only[fill1]
  mate2_start[fill2] <- mate2_start_only[fill2]
  mate1_end[fill1] <- mate1_start_only[fill1] + widths[fill1] - 1L
  mate2_end[fill2] <- mate2_start_only[fill2] + widths[fill2] - 1L

  unique(data.table(
    read_id = str_match(headers, "^([^/]+)")[, 2],
    transcript_id = str_match(headers, "/([^|;]+)\\|")[, 2],
    gene_id = str_match(headers, "\\|([^;]+)")[, 2],
    mate1_start = mate1_start,
    mate1_end = mate1_end,
    mate2_start = mate2_start,
    mate2_end = mate2_end
  ), by = "read_id")
}

load_gtf <- function(path) {
  exons <- fread(path)[
    ,
    `:=`(
      start = as.integer(start),
      end = as.integer(end),
      exon_number = as.integer(exon_number),
      exon_count = as.integer(exon_count)
    )
  ]
  setorder(exons, transcript_id, exon_number, start)

  tx_info <- split(exons, by = "transcript_id", keep.by = FALSE)
  tx_info <- lapply(tx_info, function(x) {
    widths <- x$end - x$start + 1L
    x[, `:=`(cum_start = cumsum(widths) - widths + 1L, cum_end = cumsum(widths))]
    x[, .(seqnames, start, end, strand, cum_start, cum_end)]
  })

  tx_table <- unique(exons[, .(transcript_id, gene_id, exon_count, transcript_type)])
  list(exons = exons, tx_info = tx_info, tx_table = tx_table)
}

read_bam <- function(path) {
  bam <- scanBam(path, param = ScanBamParam(what = c("qname", "rname", "pos", "cigar", "strand", "flag", "mapq")))[[1]]
  dt <- data.table(
    qname = bam$qname,
    read_id = str_match(bam$qname, "^([^/]+)")[, 2],
    rname = as.character(bam$rname),
    pos = bam$pos,
    cigar = bam$cigar,
    strand = as.character(bam$strand),
    flag = bam$flag,
    mapq = bam$mapq
  )
  dt[, `:=`(
    is_unmapped = bitwAnd(flag, 4L) != 0L,
    is_secondary = bitwAnd(flag, 256L) != 0L,
    is_supplementary = bitwAnd(flag, 2048L) != 0L,
    is_read1 = bitwAnd(flag, 64L) != 0L,
    is_read2 = bitwAnd(flag, 128L) != 0L,
    is_spliced_observed = str_detect(cigar, "N")
  )]
}

attach_truth <- function(bam_dt, truth1, truth2) {
  x <- merge(copy(bam_dt), truth1[, .(read_id, transcript_id, gene_id, mate1_start, mate1_end)], by = "read_id", all.x = TRUE, sort = FALSE)
  x <- merge(x, truth2[, .(read_id, mate2_start, mate2_end)], by = "read_id", all.x = TRUE, sort = FALSE)
  x[, `:=`(
    true_tx_start = fifelse(is_read1, mate1_start, fifelse(is_read2, mate2_start, NA_integer_)),
    true_tx_end = fifelse(is_read1, mate1_end, fifelse(is_read2, mate2_end, NA_integer_))
  )]
}

tx_to_blocks <- function(tx_start, tx_end, tx_exons) {
  if (is.null(tx_exons) || is.na(tx_start) || is.na(tx_end)) {
    return(NULL)
  }
  keep <- vector("list", nrow(tx_exons))
  k <- 0L
  for (i in seq_len(nrow(tx_exons))) {
    s <- max(tx_start, tx_exons$cum_start[i])
    e <- min(tx_end, tx_exons$cum_end[i])
    if (s > e) next
    off_s <- s - tx_exons$cum_start[i]
    off_e <- e - tx_exons$cum_start[i]
    if (tx_exons$strand[i] == "+") {
      gstart <- tx_exons$start[i] + off_s
      gend <- tx_exons$start[i] + off_e
    } else {
      gstart <- tx_exons$end[i] - off_e
      gend <- tx_exons$end[i] - off_s
    }
    k <- k + 1L
    keep[[k]] <- list(seqnames = tx_exons$seqnames[i], start = min(gstart, gend), end = max(gstart, gend))
  }
  if (k == 0L) NULL else rbindlist(keep[seq_len(k)])
}

summarize_blocks <- function(blocks) {
  if (is.null(blocks) || nrow(blocks) == 0L) {
    return(list(chr = NA_character_, start = NA_integer_, end = NA_integer_, ranges = NA_character_, n = 0L, junctions = "", key = ""))
  }
  junctions <- if (nrow(blocks) < 2L) "" else paste0(blocks$end[-nrow(blocks)], ">", blocks$start[-1], collapse = ";")
  key <- paste0(blocks$seqnames, ":", blocks$start, "-", blocks$end, collapse = ";")
  list(chr = blocks$seqnames[[1]], start = min(blocks$start), end = max(blocks$end), ranges = key, n = nrow(blocks), junctions = junctions, key = key)
}

build_observed_cache <- function(reads) {
  keys <- unique(reads[is_unmapped == FALSE, .(rname, pos, cigar)])
  if (nrow(keys) == 0L) {
    return(data.table(rname = character(), pos = integer(), cigar = character(), obs_blocks = list(), obs_chr = character(), obs_genomic_start = integer(), obs_genomic_end = integer(), obs_ranges = character(), obs_block_count = integer(), obs_junction_key = character(), obs_block_key = character()))
  }
  raw_blocks <- cigarRangesAlongReferenceSpace(cigar = keys$cigar, pos = keys$pos, ops = c("M", "=", "X"), reduce.ranges = FALSE)
  blocks <- Map(function(chr, gr) {
    if (!length(gr)) return(NULL)
    data.table(seqnames = rep(as.character(chr), length(gr)), start = start(gr), end = end(gr))
  }, keys$rname, raw_blocks)
  sums <- lapply(blocks, summarize_blocks)
  keys[, `:=`(
    obs_blocks = blocks,
    obs_chr = vapply(sums, `[[`, character(1), "chr"),
    obs_genomic_start = vapply(sums, `[[`, integer(1), "start"),
    obs_genomic_end = vapply(sums, `[[`, integer(1), "end"),
    obs_ranges = vapply(sums, `[[`, character(1), "ranges"),
    obs_block_count = vapply(sums, `[[`, integer(1), "n"),
    obs_junction_key = vapply(sums, `[[`, character(1), "junctions"),
    obs_block_key = vapply(sums, `[[`, character(1), "key")
  )]
}

build_true_cache <- function(reads, tx_info, threads) {
  keys <- unique(
    reads[!is.na(transcript_id) & !is.na(true_tx_start) & !is.na(true_tx_end), .(transcript_id, true_tx_start, true_tx_end)]
  )
  if (nrow(keys) == 0L) {
    return(data.table(transcript_id = character(), true_tx_start = integer(), true_tx_end = integer(), true_chr = character(), true_genomic_start = integer(), true_genomic_end = integer(), true_ranges = character(), true_block_count = integer(), true_junction_key = character(), true_block_key = character()))
  }
  blocks <- pmap_rows(nrow(keys), function(i) tx_to_blocks(keys$true_tx_start[i], keys$true_tx_end[i], tx_info[[keys$transcript_id[i]]]), threads)
  sums <- lapply(blocks, summarize_blocks)
  keys[, `:=`(
    true_chr = vapply(sums, `[[`, character(1), "chr"),
    true_genomic_start = vapply(sums, `[[`, integer(1), "start"),
    true_genomic_end = vapply(sums, `[[`, integer(1), "end"),
    true_ranges = vapply(sums, `[[`, character(1), "ranges"),
    true_block_count = vapply(sums, `[[`, integer(1), "n"),
    true_junction_key = vapply(sums, `[[`, character(1), "junctions"),
    true_block_key = vapply(sums, `[[`, character(1), "key")
  )]
}

blocks_within_transcript <- function(obs_blocks, tx_exons) {
  if (is.null(obs_blocks) || is.null(tx_exons) || nrow(obs_blocks) == 0L || nrow(tx_exons) == 0L) return(FALSE)
  x <- copy(obs_blocks)[, obs_id := .I]
  y <- copy(tx_exons)[, .(seqnames, start, end)]
  setkey(x, seqnames, start, end)
  setkey(y, seqnames, start, end)
  uniqueN(foverlaps(x, y, nomatch = 0L)$obs_id) == nrow(x)
}

build_overlap_cache <- function(class_tbl, obs_cache, tx_info, threads) {
  pairs <- unique(class_tbl[is_unmapped == FALSE & !is.na(transcript_id), .(transcript_id, obs_block_key)])
  if (nrow(pairs) == 0L) {
    return(data.table(transcript_id = character(), obs_block_key = character(), correct_transcript = logical()))
  }
  obs_map <- obs_cache[, .(obs_blocks = list(obs_blocks[[1]])), by = obs_block_key]
  pairs <- merge(pairs, obs_map, by = "obs_block_key", all.x = TRUE, sort = FALSE)
  pairs[, correct_transcript := unlist(pmap_rows(.N, function(i) blocks_within_transcript(obs_blocks[[i]], tx_info[[transcript_id[i]]]), threads), use.names = FALSE)]
  pairs[, .(transcript_id, obs_block_key, correct_transcript)]
}

classify_mapping <- function(tbl, gtf, threads) {
  reads <- unique(copy(tbl), by = "qname")
  obs_cache <- build_observed_cache(reads)
  true_cache <- build_true_cache(reads, gtf$tx_info, threads)

  class_tbl <- merge(reads[, .(qname, rname, pos, cigar, transcript_id, true_tx_start, true_tx_end, is_unmapped)], obs_cache[, .(rname, pos, cigar, obs_block_key, obs_junction_key, obs_chr, obs_genomic_start, obs_genomic_end, obs_ranges, obs_block_count)], by = c("rname", "pos", "cigar"), all.x = TRUE, sort = FALSE)
  class_tbl <- merge(class_tbl, true_cache, by = c("transcript_id", "true_tx_start", "true_tx_end"), all.x = TRUE, sort = FALSE)
  class_tbl[is.na(obs_block_key), `:=`(obs_block_key = "", obs_junction_key = "", obs_block_count = 0L)]
  class_tbl[is.na(true_block_key), `:=`(true_block_key = "", true_junction_key = "", true_block_count = 0L)]

  overlap_cache <- build_overlap_cache(class_tbl, obs_cache, gtf$tx_info, threads)
  class_tbl <- merge(class_tbl, overlap_cache, by = c("transcript_id", "obs_block_key"), all.x = TRUE, sort = FALSE)
  class_tbl[is.na(correct_transcript), correct_transcript := FALSE]
  class_tbl[, `:=`(
    obs_junction_count = fifelse(obs_junction_key == "", 0L, lengths(strsplit(obs_junction_key, ";", fixed = TRUE))),
    true_junction_count = fifelse(true_junction_key == "", 0L, lengths(strsplit(true_junction_key, ";", fixed = TRUE))),
    all_junctions_match = obs_junction_key == true_junction_key,
    exact_structure = obs_block_key == true_block_key,
    same_chr = !is.na(obs_chr) & !is.na(true_chr) & obs_chr == true_chr
  )]
  class_tbl[, .(qname, obs_chr, true_chr, obs_genomic_start, obs_genomic_end, true_genomic_start, true_genomic_end, obs_ranges, true_ranges, obs_block_count, true_block_count, true_junction_count, obs_junction_count, correct_transcript, all_junctions_match, exact_structure, same_chr)]
}

standard_summary <- function(tbl, sample, mapper, dataset, run) {
  total <- nrow(tbl)
  mapped <- tbl[, sum(!is_unmapped)]
  data.table(sample_id = sample, mapper = mapper, dataset = dataset, run = run, total_reads = total, mapped_reads = mapped, unmapped_reads = tbl[, sum(is_unmapped)], pct_mapped = 100 * mapped / total, pct_unmapped = 100 * sum(tbl$is_unmapped) / total, approx_unique_alignments = tbl[!is_unmapped & !is_secondary & !is_supplementary, .N], reads_with_multiple_alignments = tbl[, .N, by = read_id][N > 1L, .N])
}

ground_summary <- function(tbl, sample, mapper, dataset, run) {
  total <- nrow(tbl)
  correct <- tbl[mapping_status == "correct", .N]
  incorrect <- tbl[mapping_status == "incorrect", .N]
  precision <- if ((correct + incorrect) > 0) correct / (correct + incorrect) else NA_real_
  recall <- if (total > 0) correct / total else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA_real_
  data.table(sample_id = sample, mapper = mapper, dataset = dataset, run = run, total_reads = total, correct_reads = correct, incorrect_reads = incorrect, unmapped_reads = tbl[mapping_status == "unmapped", .N], accuracy = correct / total, precision = precision, recall = recall, f1_score = f1)
}

stratified_summary <- function(tbl, sample, mapper, dataset, run) {
  copy(tbl)[, `:=`(sample_id = sample, mapper = mapper, dataset = dataset, run = run)][, .(n_reads = .N, n_correct = sum(mapping_status == "correct"), n_incorrect = sum(mapping_status == "incorrect"), n_unmapped = sum(mapping_status == "unmapped")), by = .(sample_id, mapper, dataset, run, read_type, transcript_type, error_type)]
}

message("Parsing truth headers")
truth1 <- parse_truth_headers(opt$`truth-read1`)
truth2 <- parse_truth_headers(opt$`truth-read2`)

message("Loading GTF table")
gtf <- load_gtf(opt$`gtf-table`)

message("Reading BAM")
bam_tbl <- read_bam(opt$bam)

message("Attaching truth metadata")
full_tbl <- attach_truth(bam_tbl, truth1, truth2)

message("Classifying alignments")
classified <- classify_mapping(full_tbl, gtf, max(1L, opt$threads))
classified <- merge(merge(copy(full_tbl), gtf$tx_table, by = c("transcript_id", "gene_id"), all.x = TRUE, sort = FALSE), classified, by = "qname", all.x = TRUE, sort = FALSE)
classified[, read_type := fifelse(is_spliced_observed, "spliced", "simple")]
classified[, mapping_status := fcase(
  is_unmapped, "unmapped",
  is.na(transcript_id), "unknown_truth",
  exact_structure, "correct",
  default = "incorrect"
)]
classified[, error_type := fcase(
  is_unmapped, "unmapped",
  is.na(transcript_id), "unknown_truth",
  exact_structure, "none",
  read_type == "spliced" & correct_transcript & !all_junctions_match, "incorrect_junction",
  correct_transcript, "incorrect_position_within_region",
  !correct_transcript, "incorrect_position",
  default = "unknown"
)]
classified[, `:=`(
  sample_id = opt$sample,
  mapper = opt$mapper,
  dataset = opt$dataset,
  run = opt$run
)]

message("Writing outputs")
write_tsv(classified, opt$`out-coordinates`)
write_tsv(standard_summary(bam_tbl, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-standard`)
write_tsv(ground_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-ground`)
write_tsv(stratified_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-stratified`)
message("Done")
