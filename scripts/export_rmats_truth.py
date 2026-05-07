#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


FIELDS = (
    "dataset",
    "gene_id",
    "transcript_id",
    "fc_label",
    "cond1_mult",
    "cond2_mult",
    "truth_status",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export simulated differential-expression truth for rMATS comparison."
    )
    parser.add_argument("--design", action="append", default=[], help="Polyester transcript_design.tsv")
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def truth_rows(path: Path):
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"dataset", "transcript_id", "cond1_mult", "cond2_mult", "transcript_gene_id"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
        for row in reader:
            cond1 = float(row["cond1_mult"])
            cond2 = float(row["cond2_mult"])
            if cond1 == cond2:
                continue
            yield {
                "dataset": row["dataset"],
                "gene_id": row["transcript_gene_id"],
                "transcript_id": row["transcript_id"],
                "fc_label": row.get("fc_label", ""),
                "cond1_mult": row["cond1_mult"],
                "cond2_mult": row["cond2_mult"],
                "truth_status": "simulated_differential_isoform",
            }


def main():
    args = parse_args()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS)
        writer.writeheader()
        for design in args.design:
            for row in truth_rows(Path(design)):
                writer.writerow(row)


if __name__ == "__main__":
    main()
