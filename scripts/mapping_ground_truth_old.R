#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(Rsamtools)
  library(data.table)
  library(stringr)
  library(Biostrings)
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

write_tsv <- function(df, path) {
  fwrite(as.data.table(df), file = path, sep = "\t", quote = FALSE, na = "NA")
}

parse_truth_headers <- function(path) {
  truth_reads <- readBStringSet(path, format = "fasta")
  headers <- names(truth_reads)
  read_lengths <- width(truth_reads)
  rm(truth_reads)

  mate1_range <- str_match(headers, "mate1:([0-9]+)-([0-9]+)")
  mate2_range <- str_match(headers, "mate2:([0-9]+)-([0-9]+)")
  mate1_start_only <- as.integer(str_match(headers, "mate1Start:([0-9]+)")[, 2])
  mate2_start_only <- as.integer(str_match(headers, "mate2Start:([0-9]+)")[, 2])

  mate1_start <- as.integer(mate1_range[, 2])
  mate1_end <- as.integer(mate1_range[, 3])
  mate2_start <- as.integer(mate2_range[, 2])
  mate2_end <- as.integer(mate2_range[, 3])

  needs_mate1_end <- is.na(mate1_end) & !is.na(mate1_start_only)
  needs_mate2_end <- is.na(mate2_end) & !is.na(mate2_start_only)

  mate1_start[is.na(mate1_start)] <- mate1_start_only[is.na(mate1_start)]
  mate2_start[is.na(mate2_start)] <- mate2_start_only[is.na(mate2_start)]
  mate1_end[needs_mate1_end] <- mate1_start_only[needs_mate1_end] + read_lengths[needs_mate1_end] - 1L
  mate2_end[needs_mate2_end] <- mate2_start_only[needs_mate2_end] + read_lengths[needs_mate2_end] - 1L

  truth_dt <- data.table(header = headers)[
    ,
    `:=`(
      read_id = str_match(header, "^([^/]+)")[, 2],
      transcript_id = str_match(header, "/([^|;]+)\\|")[, 2],
      gene_id = str_match(header, "\\|([^;]+)")[, 2],
      mate1_start = mate1_start,
      mate1_end = mate1_end,
      mate2_start = mate2_start,
      mate2_end = mate2_end
    )
  ][
    ,
    .(read_id, transcript_id, gene_id, mate1_start, mate1_end, mate2_start, mate2_end)
  ]

  unique(truth_dt, by = "read_id")
}

load_gtf_tables <- function(gtf_table_file) {
  exon_dt <- fread(gtf_table_file, sep = "\t")
  setDT(exon_dt)

  exon_dt[
    ,
    `:=`(
      start = as.integer(start),
      end = as.integer(end),
      exon_number = as.integer(exon_number),
      exon_count = as.integer(exon_count)
    )
  ]
  setorder(exon_dt, transcript_id, exon_number, start)

  tx_dt <- unique(
    exon_dt[, .(transcript_id, gene_id, exon_count, transcript_type)],
    by = c("transcript_id", "gene_id", "exon_count", "transcript_type")
  )

  list(exons = exon_dt, transcripts = tx_dt)
}

build_transcript_index <- function(exons) {
  tx_exons <- split(exons, by = "transcript_id", keep.by = FALSE)

  lapply(tx_exons, function(tx_dt) {
    widths <- tx_dt$end - tx_dt$start + 1L
    cum_end <- cumsum(widths)
    cum_start <- cum_end - widths + 1L

    tx_dt <- copy(tx_dt)[
      ,
      `:=`(
        cum_start = cum_start,
        cum_end = cum_end
      )
    ]

    list(
      exons = tx_dt[, .(seqnames, start, end, strand, exon_number, cum_start, cum_end)]
    )
  })
}

bam_to_table <- function(bam_file) {
  bam <- scanBam(
    bam_file,
    param = ScanBamParam(what = c("qname", "rname", "pos", "cigar", "strand", "flag", "mapq"))
  )[[1]]

  bam_dt <- data.table(
    qname = bam$qname,
    read_id = str_match(bam$qname, "^([^/]+)")[, 2],
    rname = as.character(bam$rname),
    pos = bam$pos,
    cigar = bam$cigar,
    strand = as.character(bam$strand),
    flag = bam$flag,
    mapq = bam$mapq
  )

  bam_dt[
    ,
    `:=`(
      is_unmapped = bitwAnd(flag, 4L) != 0L,
      is_secondary = bitwAnd(flag, 256L) != 0L,
      is_supplementary = bitwAnd(flag, 2048L) != 0L,
      is_read1 = bitwAnd(flag, 64L) != 0L,
      is_read2 = bitwAnd(flag, 128L) != 0L,
      is_spliced_observed = str_detect(cigar, "N")
    )
  ]

  bam_dt
}

attach_truth <- function(bam_dt, truth1, truth2) {
  truth1_dt <- copy(truth1)[, .(read_id, transcript_id, gene_id, mate1_start, mate1_end)]
  truth2_dt <- copy(truth2)[, .(read_id, mate2_start, mate2_end)]

  merged_dt <- merge(copy(bam_dt), truth1_dt, by = "read_id", all.x = TRUE, sort = FALSE)
  merged_dt <- merge(merged_dt, truth2_dt, by = "read_id", all.x = TRUE, sort = FALSE)

  merged_dt[
    ,
    `:=`(
      true_tx_start = fifelse(is_read1, mate1_start, fifelse(is_read2, mate2_start, NA_integer_)),
      true_tx_end = fifelse(is_read1, mate1_end, fifelse(is_read2, mate2_end, NA_integer_))
    )
  ]

  merged_dt
}

tx_to_blocks <- function(tx_start, tx_end, tx_info) {
  if (is.null(tx_info) || is.na(tx_start) || is.na(tx_end)) {
    return(NULL)
  }

  tx_exons <- tx_info$exons
  keep <- vector("list", nrow(tx_exons))
  keep_idx <- 0L

  for (i in seq_len(nrow(tx_exons))) {
    s <- max(tx_start, tx_exons$cum_start[i])
    e <- min(tx_end, tx_exons$cum_end[i])
    if (s > e) {
      next
    }

    off_s <- s - tx_exons$cum_start[i]
    off_e <- e - tx_exons$cum_start[i]

    if (tx_exons$strand[i] == "+") {
      gstart <- tx_exons$start[i] + off_s
      gend <- tx_exons$start[i] + off_e
    } else {
      gstart <- tx_exons$end[i] - off_e
      gend <- tx_exons$end[i] - off_s
    }

    keep_idx <- keep_idx + 1L
    keep[[keep_idx]] <- list(
      seqnames = tx_exons$seqnames[i],
      start = min(gstart, gend),
      end = max(gstart, gend)
    )
  }

  if (keep_idx == 0L) {
    return(NULL)
  }

  rbindlist(keep[seq_len(keep_idx)])
}

summarize_blocks <- function(blocks) {
  if (is.null(blocks) || nrow(blocks) == 0L) {
    return(list(
      chr = NA_character_,
      ranges = NA_character_,
      n = 0L,
      start = NA_integer_,
      end = NA_integer_,
      junction_count = 0L,
      junction_key = "",
      block_key = ""
    ))
  }

  junctions <- if (nrow(blocks) < 2L) character() else paste0(blocks$end[-nrow(blocks)], ">", blocks$start[-1])
  range_key <- paste0(blocks$seqnames, ":", blocks$start, "-", blocks$end, collapse = ";")

  list(
    chr = blocks$seqnames[[1]],
    ranges = range_key,
    n = nrow(blocks),
    start = min(blocks$start),
    end = max(blocks$end),
    junction_count = length(junctions),
    junction_key = paste(junctions, collapse = ";"),
    block_key = range_key
  )
}

build_observed_cache <- function(unique_reads) {
  mapped_reads <- unique_reads[
    is_unmapped == FALSE,
    .(rname, pos, cigar)
  ]
  mapped_reads <- unique(mapped_reads)
  if (nrow(mapped_reads) == 0L) {
    return(data.table(
      rname = character(),
      pos = integer(),
      cigar = character(),
      obs_blocks = list(),
      obs_chr = character(),
      obs_genomic_start = integer(),
      obs_genomic_end = integer(),
      obs_ranges = character(),
      obs_block_count = integer(),
      obs_junction_count = integer(),
      obs_junction_key = character(),
      obs_block_key = character()
    ))
  }

  ranges_list <- cigarRangesAlongReferenceSpace(
    cigar = mapped_reads$cigar,
    pos = mapped_reads$pos,
    ops = c("M", "=", "X"),
    reduce.ranges = FALSE
  )

  blocks_list <- Map(function(chr, gr) {
    if (!length(gr)) {
      return(NULL)
    }

    data.table(
      seqnames = rep(as.character(chr), length(gr)),
      start = start(gr),
      end = end(gr)
    )
  }, mapped_reads$rname, ranges_list)

  summaries <- lapply(blocks_list, summarize_blocks)

  mapped_reads[
    ,
    `:=`(
      obs_blocks = blocks_list,
      obs_chr = vapply(summaries, `[[`, character(1), "chr"),
      obs_genomic_start = vapply(summaries, `[[`, integer(1), "start"),
      obs_genomic_end = vapply(summaries, `[[`, integer(1), "end"),
      obs_ranges = vapply(summaries, `[[`, character(1), "ranges"),
      obs_block_count = vapply(summaries, `[[`, integer(1), "n"),
      obs_junction_count = vapply(summaries, `[[`, integer(1), "junction_count"),
      obs_junction_key = vapply(summaries, `[[`, character(1), "junction_key"),
      obs_block_key = vapply(summaries, `[[`, character(1), "block_key")
    )
  ]
}

build_true_cache <- function(unique_reads, transcript_index) {
  truth_reads <- unique_reads[
    !is.na(transcript_id) & !is.na(true_tx_start) & !is.na(true_tx_end),
    .(transcript_id, true_tx_start, true_tx_end)
  ]
  truth_reads <- unique(truth_reads)

  if (nrow(truth_reads) == 0L) {
    return(data.table(
      transcript_id = character(),
      true_tx_start = integer(),
      true_tx_end = integer(),
      true_chr = character(),
      true_genomic_start = integer(),
      true_genomic_end = integer(),
      true_ranges = character(),
      true_block_count = integer(),
      true_junction_count = integer(),
      true_junction_key = character(),
      true_block_key = character()
    ))
  }

  blocks_list <- lapply(seq_len(nrow(truth_reads)), function(i) {
    tx_id <- truth_reads$transcript_id[[i]]
    tx_info <- transcript_index[[tx_id]]
    tx_to_blocks(truth_reads$true_tx_start[[i]], truth_reads$true_tx_end[[i]], tx_info)
  })

  summaries <- lapply(blocks_list, summarize_blocks)

  truth_reads[
    ,
    `:=`(
      true_chr = vapply(summaries, `[[`, character(1), "chr"),
      true_genomic_start = vapply(summaries, `[[`, integer(1), "start"),
      true_genomic_end = vapply(summaries, `[[`, integer(1), "end"),
      true_ranges = vapply(summaries, `[[`, character(1), "ranges"),
      true_block_count = vapply(summaries, `[[`, integer(1), "n"),
      true_junction_count = vapply(summaries, `[[`, integer(1), "junction_count"),
      true_junction_key = vapply(summaries, `[[`, character(1), "junction_key"),
      true_block_key = vapply(summaries, `[[`, character(1), "block_key")
    )
  ]
}

blocks_overlap_transcript <- function(obs_blocks, transcript_exons) {
  if (is.null(obs_blocks) || nrow(obs_blocks) == 0L || is.null(transcript_exons) || nrow(transcript_exons) == 0L) {
    return(FALSE)
  }

  obs_dt <- copy(obs_blocks)[, obs_id := .I]
  exon_dt <- copy(transcript_exons)[, .(seqnames, start, end)]
  setkey(obs_dt, seqnames, start, end)
  setkey(exon_dt, seqnames, start, end)

  overlaps <- foverlaps(obs_dt, exon_dt, nomatch = 0L)
  uniqueN(overlaps$obs_id) == nrow(obs_dt)
}

build_overlap_cache <- function(unique_reads, observed_cache, transcript_index) {
  if (nrow(observed_cache) == 0L) {
    return(data.table(
      transcript_id = character(),
      obs_block_key = character(),
      correct_transcript = logical()
    ))
  }

  obs_block_map <- observed_cache[
    ,
    .(obs_blocks = list(obs_blocks[[1]])),
    by = obs_block_key
  ]
  overlap_inputs <- unique_reads[
    is_unmapped == FALSE & !is.na(transcript_id),
    .(transcript_id, obs_block_key)
  ]
  overlap_inputs <- unique(overlap_inputs)

  if (nrow(overlap_inputs) == 0L) {
    return(data.table(
      transcript_id = character(),
      obs_block_key = character(),
      correct_transcript = logical()
    ))
  }

  overlap_inputs <- merge(overlap_inputs, obs_block_map, by = "obs_block_key", all.x = TRUE, sort = FALSE)
  overlap_inputs[
    ,
    correct_transcript := vapply(
      seq_len(.N),
      function(i) {
        tx_info <- transcript_index[[transcript_id[[i]]]]
        transcript_exons <- if (is.null(tx_info)) NULL else tx_info$exons
        blocks_overlap_transcript(obs_blocks[[i]], transcript_exons)
      },
      logical(1)
    )
  ][
    ,
    .(transcript_id, obs_block_key, correct_transcript)
  ]
}

classify_mapping <- function(tbl, exons, tx_tbl) {
  transcript_index <- build_transcript_index(exons)
  unique_reads <- unique(copy(tbl), by = "qname")
  observed_cache <- build_observed_cache(unique_reads)
  true_cache <- build_true_cache(unique_reads, transcript_index)

  class_tbl <- merge(
    unique_reads[, .(qname, rname, pos, cigar, transcript_id, true_tx_start, true_tx_end, is_unmapped)],
    observed_cache[, .(
      rname, pos, cigar, obs_block_key, obs_junction_key, obs_chr,
      obs_genomic_start, obs_genomic_end, obs_ranges, obs_block_count,
      obs_junction_count
    )],
    by = c("rname", "pos", "cigar"),
    all.x = TRUE,
    sort = FALSE
  )

  class_tbl <- merge(
    class_tbl,
    true_cache,
    by = c("transcript_id", "true_tx_start", "true_tx_end"),
    all.x = TRUE,
    sort = FALSE
  )

  class_tbl[
    is.na(obs_block_key),
    `:=`(
      obs_block_key = "",
      obs_junction_key = "",
      obs_block_count = 0L,
      obs_junction_count = 0L
    )
  ]
  class_tbl[
    is.na(true_block_key),
    `:=`(
      true_block_key = "",
      true_junction_key = "",
      true_block_count = 0L,
      true_junction_count = 0L
    )
  ]

  overlap_cache <- build_overlap_cache(class_tbl, observed_cache, transcript_index)
  class_tbl <- merge(class_tbl, overlap_cache, by = c("transcript_id", "obs_block_key"), all.x = TRUE, sort = FALSE)
  class_tbl[is.na(correct_transcript), correct_transcript := FALSE]

  class_tbl[
    ,
    `:=`(
      exact_structure = obs_block_key == true_block_key,
      all_junctions_match = obs_junction_key == true_junction_key,
      same_chr = !is.na(obs_chr) & !is.na(true_chr) & obs_chr == true_chr
    )
  ]

  class_tbl <- class_tbl[
    ,
    .(
      qname,
      obs_chr,
      true_chr,
      obs_genomic_start,
      obs_genomic_end,
      true_genomic_start,
      true_genomic_end,
      obs_ranges,
      true_ranges,
      obs_block_count,
      true_block_count,
      true_junction_count,
      obs_junction_count,
      correct_transcript,
      all_junctions_match,
      exact_structure,
      same_chr
    )
  ]

  classified_dt <- merge(copy(tbl), tx_tbl, by = c("transcript_id", "gene_id"), all.x = TRUE, sort = FALSE)
  classified_dt <- merge(classified_dt, class_tbl, by = "qname", all.x = TRUE, sort = FALSE)

  classified_dt[, read_type := fifelse(is_spliced_observed, "spliced", "simple")]
  classified_dt[
    ,
    mapping_status := fcase(
      is_unmapped, "unmapped",
      is.na(transcript_id), "unknown_truth",
      exact_structure, "correct",
      default = "incorrect"
    )
  ]
  classified_dt[
    ,
    error_type := fcase(
      is_unmapped, "unmapped",
      is.na(transcript_id), "unknown_truth",
      exact_structure, "none",
      read_type == "spliced" & correct_transcript & !all_junctions_match, "incorrect_junction",
      correct_transcript, "incorrect_position_within_region",
      !correct_transcript, "incorrect_position",
      default = "unknown"
    )
  ]

  classified_dt
}

standard_summary <- function(tbl, sample, mapper, dataset, run) {
  total <- nrow(tbl)
  mapped <- tbl[, sum(!is_unmapped)]
  unmapped <- tbl[, sum(is_unmapped)]
  unique_aln <- tbl[!is_unmapped & !is_secondary & !is_supplementary, .N]
  multi <- tbl[, .N, by = read_id][N > 1L, .N]

  data.table(sample_id = sample,
    mapper = mapper, dataset = dataset, run = run,
    total_reads = total, mapped_reads = mapped, unmapped_reads = unmapped,
    pct_mapped = 100 * mapped / total, pct_unmapped = 100 * unmapped / total,
    approx_unique_alignments = unique_aln, reads_with_multiple_alignments = multi
  )
}

ground_summary <- function(tbl, sample, mapper, dataset, run) {
  total <- nrow(tbl)
  correct <- tbl[mapping_status == "correct", .N]
  incorrect <- tbl[mapping_status == "incorrect", .N]
  unmapped <- tbl[mapping_status == "unmapped", .N]

  precision <- if ((correct + incorrect) > 0) correct / (correct + incorrect) else NA_real_
  recall <- if (total > 0) correct / total else NA_real_
  f1 <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }

  data.table(sample_id = sample, mapper = mapper,
    dataset = dataset, run = run, total_reads = total, correct_reads = correct,
    incorrect_reads = incorrect, unmapped_reads = unmapped,
    accuracy = correct / total, precision = precision,
    recall = recall, f1_score = f1)
}

stratified_summary <- function(tbl, sample, mapper, dataset, run) {
  summary_dt <- copy(tbl)
  summary_dt[
    ,
    `:=`(
      sample_id = sample, mapper = mapper,
      dataset = dataset, run = run
    )]

  summary_dt[
    ,
    .(n_reads = .N,
      n_correct = sum(mapping_status == "correct"),
      n_incorrect = sum(mapping_status == "incorrect"),
      n_unmapped = sum(mapping_status == "unmapped")
    ),
    by = .(sample_id, mapper, dataset, run, read_type, transcript_type, error_type)
  ]
}

message("Parsing truth headers")
truth1 <- parse_truth_headers(opt$`truth-read1`)
truth2 <- parse_truth_headers(opt$`truth-read2`)

message("Loading GTF table")
gtf <- load_gtf_tables(opt$`gtf-table`)

message("Reading BAM")
bam_tbl <- bam_to_table(opt$bam)

message("Attaching truth metadata")
full_tbl <- attach_truth(bam_tbl, truth1, truth2)

message("Classifying alignments")
classified <- classify_mapping(full_tbl, gtf$exons, gtf$transcripts)
classified[, `:=`(sample_id = opt$sample,mapper = opt$mapper,dataset = opt$dataset,run = opt$run)]

message("Writing outputs")
write_tsv(classified, opt$`out-coordinates`)
write_tsv(standard_summary(bam_tbl, opt$sample, opt$mapper, opt$dataset, opt$run),opt$`out-standard`)
write_tsv(ground_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run),opt$`out-ground`)
write_tsv(stratified_summary(classified, opt$sample, opt$mapper, opt$dataset, opt$run), opt$`out-stratified`)
message("Done")
