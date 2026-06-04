#!/usr/bin/env python3
"""Parse STAR Log.final.out or HISAT2 stderr into a one-row mapping-QC TSV.

The two mappers produce very different stat files, so we normalize them to a
common schema (uniquely-mapped / multi-mapped / unmapped counts and percentages)
plus mapper-specific columns. Missing values are emitted as empty strings so
the downstream merge in pandas treats them as NaN.
"""

import argparse
import csv
import re
from pathlib import Path


COMMON_FIELDS = [
    "sample_id", "mapper", "dataset", "run",
    "total_input_reads",
    "uniquely_mapped_reads", "pct_uniquely_mapped",
    "multi_mapped_reads", "pct_multi_mapped",
    "unmapped_reads", "pct_unmapped",
]

STAR_FIELDS = [
    "avg_input_read_length",
    "avg_mapped_length",
    "num_splices_total",
    "num_splices_annotated",
    "num_splices_gt_ag",
    "mismatch_rate_pct",
    "deletion_rate_pct",
    "insertion_rate_pct",
    "too_many_loci_reads",
    "pct_too_many_loci",
]

HISAT2_FIELDS = [
    "pct_concordant_unique",
    "pct_concordant_multi",
    "pct_discordant_unique",
    "overall_alignment_rate",
]

ALL_FIELDS = COMMON_FIELDS + STAR_FIELDS + HISAT2_FIELDS


def parse_args():
    arg_parser = argparse.ArgumentParser(description=__doc__)
    arg_parser.add_argument("--mapper", required=True, choices=["star", "hisat2"])
    arg_parser.add_argument("--log", required=True, help="STAR Log.final.out or HISAT2 stderr log path.")
    arg_parser.add_argument("--sample", required=True)
    arg_parser.add_argument("--dataset", required=True)
    arg_parser.add_argument("--run", required=True)
    arg_parser.add_argument("--out", required=True)
    return arg_parser.parse_args()


def to_float(value):
    if value is None:
        return None
    value = value.strip().rstrip("%")
    try:
        return float(value)
    except ValueError:
        return None


def to_int(value):
    float_val = to_float(value)
    return int(float_val) if float_val is not None else None


STAR_LABELS = {
    "Number of input reads": ("total_input_reads", to_int),
    "Average input read length": ("avg_input_read_length", to_float),
    "Uniquely mapped reads number": ("uniquely_mapped_reads", to_int),
    "Uniquely mapped reads %": ("pct_uniquely_mapped", to_float),
    "Average mapped length": ("avg_mapped_length", to_float),
    "Number of splices: Total": ("num_splices_total", to_int),
    "Number of splices: Annotated (sjdb)": ("num_splices_annotated", to_int),
    "Number of splices: GT/AG": ("num_splices_gt_ag", to_int),
    "Mismatch rate per base, %": ("mismatch_rate_pct", to_float),
    "Deletion rate per base": ("deletion_rate_pct", to_float),
    "Insertion rate per base": ("insertion_rate_pct", to_float),
    "Number of reads mapped to multiple loci": ("multi_mapped_reads", to_int),
    "% of reads mapped to multiple loci": ("pct_multi_mapped", to_float),
    "Number of reads mapped to too many loci": ("too_many_loci_reads", to_int),
    "% of reads mapped to too many loci": ("pct_too_many_loci", to_float),
    "% of reads unmapped: too many mismatches": ("_unmap_mismatch", to_float),
    "% of reads unmapped: too short": ("_unmap_short", to_float),
    "% of reads unmapped: other": ("_unmap_other", to_float),
}


def parse_star_log(path):
    out = {}
    text = Path(path).read_text(errors="ignore")
    for line in text.splitlines():
        if "|" not in line:
            continue
        label, value = line.split("|", 1)
        label = label.strip()
        value = value.strip()
        if label in STAR_LABELS:
            key, conv = STAR_LABELS[label]
            out[key] = conv(value)

    unmap_keys = ["_unmap_mismatch", "_unmap_short", "_unmap_other"]
    if any(unmap_key in out for unmap_key in unmap_keys):
        out["pct_unmapped"] = sum(out.get(unmap_key, 0) or 0 for unmap_key in unmap_keys)
        total = out.get("total_input_reads")
        if total and out["pct_unmapped"] is not None:
            out["unmapped_reads"] = int(round(total * out["pct_unmapped"] / 100))
    for unmap_key in unmap_keys:
        out.pop(unmap_key, None)

    if out.get("pct_too_many_loci") is not None:
        out["pct_multi_mapped"] = (out.get("pct_multi_mapped") or 0) + out["pct_too_many_loci"]
        if out.get("too_many_loci_reads") is not None:
            out["multi_mapped_reads"] = (out.get("multi_mapped_reads") or 0) + out["too_many_loci_reads"]

    return out


HISAT2_PATTERNS = [
    (re.compile(r"^\s*(\d+)\s+reads;\s+of these:"), "total_input_reads", to_int),
    (re.compile(r"^\s*\d+\s+\(([\d.]+)%\)\s+aligned concordantly exactly 1 time"),
        "pct_concordant_unique", to_float),
    (re.compile(r"^\s*\d+\s+\(([\d.]+)%\)\s+aligned concordantly >1 times"),
        "pct_concordant_multi", to_float),
    (re.compile(r"^\s*\d+\s+\(([\d.]+)%\)\s+aligned discordantly 1 time"),
        "pct_discordant_unique", to_float),
    (re.compile(r"^([\d.]+)%\s+overall alignment rate"),
        "overall_alignment_rate", to_float),
]


def parse_hisat2_log(path):
    out = {}
    text = Path(path).read_text(errors="ignore")
    for line in text.splitlines():
        for pattern, key, conv in HISAT2_PATTERNS:
            match = pattern.search(line)
            if match:
                out[key] = conv(match.group(1))
                break

    total = out.get("total_input_reads")
    if out.get("overall_alignment_rate") is not None:
        out["pct_unmapped"] = 100 - out["overall_alignment_rate"]
        if total:
            out["unmapped_reads"] = int(round(total * out["pct_unmapped"] / 100))
    if out.get("pct_concordant_unique") is not None:
        out["pct_uniquely_mapped"] = out["pct_concordant_unique"]
        if total:
            out["uniquely_mapped_reads"] = int(round(total * out["pct_uniquely_mapped"] / 100))
    if out.get("pct_concordant_multi") is not None:
        # Include discordant pairs so multi-mapped covers all non-uniquely
        # mapped reads, similar to STAR's multi-mapped category.
        disc = out.get("pct_discordant_unique") or 0
        out["pct_multi_mapped"] = out["pct_concordant_multi"] + disc
        if total:
            out["multi_mapped_reads"] = int(round(total * out["pct_multi_mapped"] / 100))
    return out


def main():
    args = parse_args()
    log_path = Path(args.log)
    parsed = {}
    if log_path.exists():
        parser = parse_star_log if args.mapper == "star" else parse_hisat2_log
        parsed = parser(log_path)

    row = {f: "" for f in ALL_FIELDS}
    row.update({
        "sample_id": args.sample,
        "mapper": args.mapper,
        "dataset": args.dataset,
        "run": args.run,
    })
    for k, v in parsed.items():
        if k in row and v is not None:
            row[k] = v

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=ALL_FIELDS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
