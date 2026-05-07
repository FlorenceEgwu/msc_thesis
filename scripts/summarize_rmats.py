#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Optional


EVENT_TYPES = ("SE", "MXE", "A5SS", "A3SS", "RI")
RMATS_SUFFIXES = (
    "MATS.JC.txt",
    "MATS.JCEC.txt",
)
SUMMARY_FIELDS = (
    "dataset",
    "mapper",
    "param_group",
    "event_type",
    "total_events",
    "significant_events",
    "significant_genes",
    "truth_supported_events",
    "truth_supported_genes",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize significant rMATS events for one dataset/mapper/parameter case."
    )
    parser.add_argument("--rmats-out", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--mapper", required=True)
    parser.add_argument("--param-group", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-abs-inc-diff", type=float, default=0.0)
    parser.add_argument("--significant-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--truth-table", default="")
    return parser.parse_args()


def rmats_file(rmats_out: Path, event_type: str) -> Optional[Path]:
    for suffix in RMATS_SUFFIXES:
        candidate = rmats_out / f"{event_type}.{suffix}"
        if candidate.exists():
            return candidate
    return None


def as_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def read_truth_genes(path: str) -> set:
    if not path:
        return set()
    truth_path = Path(path)
    if not truth_path.exists():
        return set()
    with truth_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return set()
        for column in ("gene_id", "GeneID", "geneSymbol"):
            if column in reader.fieldnames:
                return {
                    row[column]
                    for row in reader
                    if row.get(column) not in ("", None)
                }
    return set()


def summarize_event(args, event_type: str, truth_genes: set):
    source = rmats_file(Path(args.rmats_out), event_type)
    if source is None:
        return [], [], {
            "dataset": args.dataset,
            "mapper": args.mapper,
            "param_group": args.param_group,
            "event_type": event_type,
            "total_events": 0,
            "significant_events": 0,
            "significant_genes": 0,
            "truth_supported_events": 0,
            "truth_supported_genes": 0,
        }

    with source.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            fieldnames = []
            rows = []
        else:
            fieldnames = list(reader.fieldnames)
            rows = list(reader)

    if rows and ("FDR" not in fieldnames or "IncLevelDifference" not in fieldnames):
        raise ValueError(f"{source} is missing FDR or IncLevelDifference columns")

    significant = []
    significant_genes = set()
    truth_supported_genes = set()
    truth_supported_events = 0
    for row in rows:
        fdr = as_float(row.get("FDR", "nan"))
        inc_diff = abs(as_float(row.get("IncLevelDifference", "nan")))
        if fdr <= args.fdr and inc_diff >= args.min_abs_inc_diff:
            out_row = {
                "dataset": args.dataset,
                "mapper": args.mapper,
                "param_group": args.param_group,
                "event_type": event_type,
            }
            out_row.update(row)
            significant.append(out_row)
            gene_id = row.get("GeneID", "")
            if gene_id:
                significant_genes.add(gene_id)
                if gene_id in truth_genes:
                    truth_supported_events += 1
                    truth_supported_genes.add(gene_id)

    summary = {
        "dataset": args.dataset,
        "mapper": args.mapper,
        "param_group": args.param_group,
        "event_type": event_type,
        "total_events": len(rows),
        "significant_events": len(significant),
        "significant_genes": len(significant_genes),
        "truth_supported_events": truth_supported_events,
        "truth_supported_genes": len(truth_supported_genes),
    }
    return significant, fieldnames, summary


def write_significant(path: Path, rows: list, rmats_fields: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["dataset", "mapper", "param_group", "event_type"] + [
        field for field in rmats_fields
        if field not in {"dataset", "mapper", "param_group", "event_type"}
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def write_summary(path: Path, rows: list):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    truth_genes = read_truth_genes(args.truth_table)
    significant_rows = []
    significant_fields = []
    summary_rows = []

    for event_type in EVENT_TYPES:
        rows, fields, summary = summarize_event(args, event_type, truth_genes)
        significant_rows.extend(rows)
        for field in fields:
            if field not in significant_fields:
                significant_fields.append(field)
        summary_rows.append(summary)

    write_significant(Path(args.significant_out), significant_rows, significant_fields)
    write_summary(Path(args.summary_out), summary_rows)


if __name__ == "__main__":
    main()
