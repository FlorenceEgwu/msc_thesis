#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse); library(GenomicAlignments); library(Rsamtools)
  library(data.table); library(stringr); library(Biostrings)
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
  make_option("--out-coordinates", type = "character"),
  make_option("--out-standard", type = "character"),
  make_option("--out-ground", type = "character"),
  make_option("--out-stratified", type = "character")
)))

parse_truth_headers <- function(path) {
  seqlens <- fasta.seqlengths(path)
  headers <- names(seqlens)
  m1 <- str_match(headers, "mate1:([0-9]+)-([0-9]+)")
  m2 <- str_match(headers, "mate2:([0-9]+)-([0-9]+)")

  # Polyester quirk: mate2_start in header is shifted by +1 vs actual read start.
  # Subtract 1 so both mates represent 1-based inclusive 100-bp ranges.
  unique(data.table(
    read_id = sub("/.*$", "", headers),
    transcript_id = str_match(headers, "/([^|;]+)\\|")[, 2],
    gene_id = str_match(headers, "\\|([^;]+)")[, 2],
    mate1_start = as.integer(m1[, 2]),
    mate1_end = as.integer(m1[, 3]),
    mate2_start = as.integer(m2[, 2]) - 1L,
    mate2_end = as.integer(m2[, 3])
  ), by = "read_id")
}

load_flat_exons <- function(gtf_table_file) {
  exon_dt <- fread(gtf_table_file, sep = "\t")
  exon_dt[, `:=`(
    start = as.integer(start), end = as.integer(end),
    exon_number = as.integer(exon_number),
    exon_count = as.integer(exon_count),
    seqnames = as.character(seqnames),
    strand = as.character(strand)
  )]
  setorder(exon_dt, transcript_id, exon_number, start)
  exon_dt[, w := end - start + 1L]
  exon_dt[, cum_end := cumsum(w), by = transcript_id]
  exon_dt[, cum_start := cum_end - w + 1L][, w := NULL]

  as_cols <- intersect(c("as_event_type"), names(exon_dt))
  tx_dt <- unique(exon_dt[, c(
    "transcript_id", "gene_id", "exon_count", "transcript_type", as_cols
  ), with = FALSE])
  if (!"as_event_type" %in% names(tx_dt)) {
    tx_dt[, as_event_type := "unknown"]
  }
  list(exons = exon_dt, transcripts = tx_dt)
}

bam_to_table <- function(bam_file) {
  bam <- scanBam(bam_file,
    param = ScanBamParam(what = c("qname", "rname", "pos", "cigar", "flag")))[[1]]

  data.table(
    qname = bam$qname,
    read_id = sub("/.*$", "", bam$qname),
    rname = as.character(bam$rname),
    pos = bam$pos,
    cigar = bam$cigar,
    flag = bam$flag,
    is_unmapped = bitwAnd(bam$flag, 4L)    != 0L,
    is_secondary = bitwAnd(bam$flag, 256L)  != 0L,
    is_supplementary = bitwAnd(bam$flag, 2048L) != 0L,
    is_read1 = bitwAnd(bam$flag, 64L)   != 0L,
    is_read2 = bitwAnd(bam$flag, 128L)  != 0L,
    is_spliced_observed = grepl("N", bam$cigar, fixed = TRUE)
  )
}

# Compute genomic blocks from CIGAR for each unique (rname, pos, cigar) triple.
build_observed_cache <- function(bam_tbl) {
  mapped <- unique(bam_tbl[is_unmapped == FALSE, .(rname, pos, cigar)])
  mapped[, pair_id := .I]

  rng <- cigarRangesAlongReferenceSpace(
    cigar = mapped$cigar, pos = mapped$pos,
    ops = c("M", "=", "X"), reduce.ranges = FALSE)
  widths <- elementNROWS(rng)

  flat <- data.table(
    pair_id = rep(mapped$pair_id, widths),
    seqnames = rep(mapped$rname,   widths),
    start = as.integer(start(unlist(rng, use.names = FALSE))),
    end = as.integer(end  (unlist(rng, use.names = FALSE))))
  setorder(flat, pair_id, seqnames, start)

  by_pair <- flat[, {
    rk <- paste0(seqnames, ":", start, "-", end, collapse = ";")
    jk <- if (.N > 1L) paste0(head(end, -1L), ">", tail(start, -1L), collapse = ";") else ""
    .(obs_chr = seqnames[1L],
      obs_genomic_start = min(start),   obs_genomic_end  = max(end),
      obs_block_count = .N,           obs_ranges       = rk,
      obs_junction_count = as.integer(.N - 1L),
      obs_junction_key = jk,           obs_block_key    = rk)
  }, by = pair_id]

  by_key <- mapped[by_pair, on = "pair_id"][, pair_id := NULL]

  reps <- by_pair[, .(pair_id = pair_id[1L]), by = obs_block_key]
  blocks_by_key <- merge(reps, flat, by = "pair_id", sort = FALSE)[,
    .(obs_block_key, seqnames, start, end)]
  list(by_key = by_key, blocks_by_key = blocks_by_key)
}

# Convert each unique (transcript_id, tx_start, tx_end) to genomic blocks via foverlaps
# (vectorized; single pass across the whole truth set).
build_true_cache <- function(truth_inputs, flat_exons) {
  truth_reads <- unique(truth_inputs[!is.na(transcript_id) & !is.na(tx_start) & !is.na(tx_end)])
  truth_reads[, pair_id := .I]

  tq  <- truth_reads[, .(transcript_id, start = tx_start, end = tx_end, pair_id)]
  fex <- flat_exons[, .(transcript_id, seqnames, strand,
                        gstart = start, gend = end, cum_start, cum_end)]
  setnames(fex, c("cum_start", "cum_end"), c("start", "end"))
  setkey(tq,  transcript_id, start, end)
  setkey(fex, transcript_id, start, end)

  hits <- foverlaps(tq, fex, type = "any", nomatch = NULL)
  setnames(hits, c("start", "end", "i.start", "i.end"),
                 c("ex_cs", "ex_ce", "tx_s", "tx_e"))

  hits[, `:=`(
    off_s = pmax(tx_s, ex_cs) - ex_cs,
    off_e = pmin(tx_e, ex_ce) - ex_cs
  )]
  hits[, `:=`(
    blk_lo = pmin(fifelse(strand == "+", gstart + off_s, gend - off_e),
                  fifelse(strand == "+", gstart + off_e, gend - off_s)),
    blk_hi = pmax(fifelse(strand == "+", gstart + off_s, gend - off_e),
                  fifelse(strand == "+", gstart + off_e, gend - off_s))
  )]
  setorder(hits, pair_id, seqnames, blk_lo)

  summaries <- hits[, {
    rk <- paste0(seqnames, ":", blk_lo, "-", blk_hi, collapse = ";")
    jk <- if (.N > 1L) paste0(head(blk_hi, -1L), ">", tail(blk_lo, -1L), collapse = ";") else ""
    .(true_chr = seqnames[1L],
      true_genomic_start = min(blk_lo), true_genomic_end = max(blk_hi),
      true_block_count = .N,          
      true_ranges = rk,
      true_junction_count = as.integer(.N - 1L),
      true_junction_key = jk,          
      true_block_key = rk)
  }, by = pair_id]

  truth_reads[summaries, on = "pair_id"][, pair_id := NULL]
}

# "correct_transcript": every observed block overlaps at least one exon of the candidate transcript.
build_overlap_cache <- function(observed_cache, class_tbl, flat_exons) {
  pairs <- unique(class_tbl[is_unmapped == FALSE & !is.na(transcript_id) & obs_block_key != "",
                            .(transcript_id, obs_block_key)])

  blocks_numbered <- copy(observed_cache$blocks_by_key)
  blocks_numbered[, obs_id := seq_len(.N), by = obs_block_key]
  exploded <- merge(pairs, blocks_numbered, by = "obs_block_key", allow.cartesian = TRUE)
  fex <- flat_exons[, .(transcript_id, seqnames, start, end)]
  setkey(fex, transcript_id, seqnames, start, end)
  setkey(exploded, transcript_id, seqnames, start, end)

  hits <- foverlaps(exploded, fex, type = "any", nomatch = NULL)
  covered <- unique(hits[, .(transcript_id, obs_block_key, obs_id)])
  covered_counts <- covered[, .N, by = .(transcript_id, obs_block_key)]
  total_counts <- blocks_numbered[, .N, by = obs_block_key]

  res <- copy(pairs)
  res[total_counts, on = "obs_block_key", n_total   := N]
  res[covered_counts, on = c("transcript_id", "obs_block_key"), n_covered := N]
  res[is.na(n_covered), n_covered := 0L]
  res[, .(transcript_id, obs_block_key, correct_transcript = n_covered == n_total)]
}

classify_mapping <- function(tbl, flat_exons, tx_tbl) {
  tbl <- copy(tbl)[, .aln_id := .I]
  observed_cache <- build_observed_cache(tbl)

  # Truth: compute for BOTH mates (Polyester randomizes orientation ~50/50).
  true_cache <- build_true_cache(rbind(
    tbl[!is.na(mate1_start), .(transcript_id, tx_start = mate1_start, tx_end = mate1_end)],
    tbl[!is.na(mate2_start), .(transcript_id, tx_start = mate2_start, tx_end = mate2_end)]
  ), flat_exons)

  class_tbl <- merge(
    tbl[, .(.aln_id, rname, pos, cigar, transcript_id,
            mate1_start, mate1_end, mate2_start, mate2_end,
            is_unmapped, is_spliced_observed)],
    observed_cache$by_key,
    by = c("rname", "pos", "cigar"), all.x = TRUE, sort = FALSE)

  true_cols <- c("true_chr","true_genomic_start","true_genomic_end",
                 "true_block_count","true_ranges","true_junction_count",
                 "true_junction_key","true_block_key")
  rename_for <- function(tc, prefix) {
    out <- copy(tc)
    setnames(out, true_cols, paste0(prefix, "_", sub("^true_", "", true_cols)))
    out
  }

  class_tbl <- merge(class_tbl, rename_for(true_cache, "m1"),
    by.x = c("transcript_id","mate1_start","mate1_end"),
    by.y = c("transcript_id","tx_start","tx_end"),
    all.x = TRUE, sort = FALSE)
  class_tbl <- merge(class_tbl, rename_for(true_cache, "m2"),
    by.x = c("transcript_id","mate2_start","mate2_end"),
    by.y = c("transcript_id","tx_start","tx_end"),
    all.x = TRUE, sort = FALSE)

  # Fill NAs with empty strings / zeros so equality tests are well-defined.
  for (col in c("obs_block_key","obs_junction_key")) set(class_tbl, which(is.na(class_tbl[[col]])), col, "")
  for (col in c("obs_block_count","obs_junction_count")) set(class_tbl, which(is.na(class_tbl[[col]])), col, 0L)
  for (p in c("m1","m2")) {
    for (col in paste0(p, c("_block_key","_junction_key"))) set(class_tbl, which(is.na(class_tbl[[col]])), col, "")
    for (col in paste0(p, c("_block_count","_junction_count"))) set(class_tbl, which(is.na(class_tbl[[col]])), col, 0L)
  }

  overlap_cache <- build_overlap_cache(observed_cache, class_tbl, flat_exons)
  class_tbl <- merge(class_tbl, overlap_cache,
    by = c("transcript_id","obs_block_key"), all.x = TRUE, sort = FALSE)
  class_tbl[is.na(correct_transcript), correct_transcript := FALSE]

  # Exact structure: obs matches EITHER mate's block layout. Pick matching mate's
  # truth for output; if neither matches, prefer m1 (else m2) for reporting.
  class_tbl[, `:=`(
    exact_m1 = obs_block_key != "" & m1_block_key != "" & obs_block_key == m1_block_key,
    exact_m2 = obs_block_key != "" & m2_block_key != "" & obs_block_key == m2_block_key
  )]
  class_tbl[, exact_structure := exact_m1 | exact_m2]
  class_tbl[, use_m2 := (!exact_m1 & exact_m2) |
                        (!exact_m1 & !exact_m2 & m1_block_key == "" & m2_block_key != "")]
  for (field in c("chr","genomic_start","genomic_end","ranges","block_count",
                  "block_key","junction_count","junction_key")) {
    m1c <- paste0("m1_", field); m2c <- paste0("m2_", field); tgt <- paste0("true_", field)
    class_tbl[, (tgt) := fifelse(use_m2, get(m2c), get(m1c))]
  }
  class_tbl[, `:=`(
    all_junctions_match = obs_junction_key == true_junction_key,
    same_chr            = !is.na(obs_chr) & !is.na(true_chr) & obs_chr == true_chr
  )]

  classified_dt <- merge(tbl, tx_tbl, by = c("transcript_id","gene_id"),
                         all.x = TRUE, sort = FALSE)
  classified_dt <- merge(classified_dt,
    class_tbl[, .(.aln_id, obs_chr, true_chr,
                  obs_genomic_start, obs_genomic_end,
                  true_genomic_start, true_genomic_end,
                  obs_ranges, true_ranges,
                  obs_block_count, true_block_count,
                  true_junction_count, obs_junction_count,
                  correct_transcript, all_junctions_match,
                  exact_structure, same_chr)],
    by = ".aln_id", all.x = TRUE, sort = FALSE)
  classified_dt[, .aln_id := NULL]

  classified_dt[, read_type := fifelse(is_spliced_observed, "spliced", "simple")]
  classified_dt[, mapping_status := fcase(
    is_unmapped, "unmapped",
    is.na(transcript_id), "unknown_truth",
    exact_structure, "correct",
    default = "incorrect"
  )]
  classified_dt[, error_type := fcase(
    is_unmapped, "unmapped",
    is.na(transcript_id), "unknown_truth",
    exact_structure, "none",
    correct_transcript & (is_spliced_observed | true_junction_count > 0L) & !all_junctions_match,
      "incorrect_junction",
    correct_transcript, "incorrect_position_within_region",
    !correct_transcript, "incorrect_position",
    default = "unknown"
  )]
  classified_dt
}

standard_summary <- function(tbl, sample, mapper, dataset, run) {
  primary <- tbl[is_secondary == FALSE & is_supplementary == FALSE]
  n_tot <- primary[, uniqueN(paste(qname, is_read1))]
  n_map <- primary[is_unmapped == FALSE, uniqueN(paste(qname, is_read1))]
  n_unmap <- primary[is_unmapped == TRUE,  uniqueN(paste(qname, is_read1))]
  n_multi <- tbl[is_unmapped == FALSE, .N, by = .(qname, is_read1)][N > 1L, .N]
  
  data.table(sample_id = sample, mapper = mapper, dataset = dataset, run = run,
             total_reads = n_tot,
             mapped_reads = n_map,
             unmapped_reads = n_unmap,
             pct_mapped = if (n_tot > 0) 100 * n_map / n_tot else NA_real_,
             reads_with_multiple_alignments = n_multi)
}

ground_summary <- function(tbl, sample, mapper, dataset, run) {
  p <- tbl[is_secondary == FALSE & is_supplementary == FALSE]
  n_tot <- nrow(p)
  n_cor <- p[mapping_status == "correct", .N]
  n_incor <- p[mapping_status == "incorrect", .N]
  n_unmap <- p[mapping_status == "unmapped", .N]
  prec <- if (n_cor + n_incor > 0) n_cor / (n_cor + n_incor) else NA_real_
  rec <- if (n_tot > 0) n_cor / n_tot  else NA_real_
  f1 <- if (!is.na(prec) && !is.na(rec) && prec + rec > 0) 2*prec*rec/(prec+rec) else NA_real_
  
  data.table(sample_id = sample, mapper = mapper, dataset = dataset, run = run,
             total_reads = n_tot,
             correct_reads = n_cor,
             incorrect_reads = n_incor,
             unmapped_reads = n_unmap,
             accuracy  = if (n_tot > 0) n_cor / n_tot else NA_real_,
             precision = prec, recall = rec, f1_score = f1)
}

stratified_summary <- function(tbl, sample, mapper, dataset, run) {
  p <- tbl[is_secondary == FALSE & is_supplementary == FALSE]
  if (!"as_event_type" %in% names(p)) {
    p[, as_event_type := "unknown"]
  }
  p[, `:=`(sample_id = sample, mapper = mapper, dataset = dataset, run = run)]
  p[, .(n_reads = .N,
        n_correct = sum(mapping_status == "correct"),
        n_incorrect = sum(mapping_status == "incorrect"),
        n_unmapped = sum(mapping_status == "unmapped")),
    by = .(sample_id, mapper, dataset, run, read_type, transcript_type, as_event_type, error_type)]
}

message("Parsing truth headers"); truth1 <- parse_truth_headers(opt$`truth-read1`)
                                  truth2 <- parse_truth_headers(opt$`truth-read2`)
message("Loading GTF table"); gtf <- load_flat_exons(opt$`gtf-table`)
message("Reading BAM"); bam_tbl <- bam_to_table(opt$bam)

message("Attaching truth metadata")
full_tbl <- merge(bam_tbl,
  truth1[, .(read_id, transcript_id, gene_id, mate1_start, mate1_end)],
  by = "read_id", all.x = TRUE, sort = FALSE)
full_tbl <- merge(full_tbl, truth2[, .(read_id, mate2_start, mate2_end)],
  by = "read_id", all.x = TRUE, sort = FALSE)

message("Classifying alignments")
classified <- classify_mapping(full_tbl, gtf$exons, gtf$transcripts)
classified[, `:=`(sample_id = opt$sample, mapper = opt$mapper,
                  dataset = opt$dataset, run = opt$run)]

message("Writing outputs")
fwrite(classified, opt$`out-coordinates`, sep = "\t", na = "NA", quote = FALSE)
fwrite(standard_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-standard`, sep = "\t", na = "NA", quote = FALSE)
fwrite(ground_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-ground`, sep = "\t", na = "NA", quote = FALSE)
fwrite(stratified_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-stratified`,  sep = "\t", na = "NA", quote = FALSE)
message("Done")
