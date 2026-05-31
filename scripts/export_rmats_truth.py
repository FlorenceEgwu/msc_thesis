#!/usr/bin/env python3
"""Export simulated differential-expression truth for rMATS comparison.

Produces two files:

- --out             per-transcript truth: one row per transcript with cond1 != cond2
- --out-genes       per-gene truth: aggregates transcript directions into a class
                    ('none', 'single_shift', 'isoform_switch', 'co_directional_shift')
                    that lets `summarize_rmats.py` distinguish true isoform-switch
                    signal from gene-level co-directional shifts and pure null genes
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


TX_FIELDS = (
    "dataset",
    "gene_id",
    "transcript_id",
    "fc_label",
    "cond1_mult",
    "cond2_mult",
    "transcript_direction",
)

GENE_FIELDS = (
    "dataset",
    "gene_id",
    "n_transcripts",
    "n_transcripts_diff",
    "n_up_in_cond1",
    "n_up_in_cond2",
    "gene_truth_class",
    "up_in_cond1_transcripts",
    "up_in_cond2_transcripts",
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", action="append", default=[],
                        help="Polyester transcript_design.tsv (repeatable, one per dataset).")
    parser.add_argument("--out", required=True, help="Per-transcript truth TSV.")
    parser.add_argument("--out-genes", required=True, help="Per-gene truth TSV.")
    return parser.parse_args()


def read_design(path: Path) -> pd.DataFrame:
    required = {"dataset", "transcript_id", "cond1_mult", "cond2_mult", "transcript_gene_id"}
    df = pd.read_csv(path, sep="\t", dtype={"cond1_mult": float, "cond2_mult": float})
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {', '.join(sorted(missing))}")
    df = df.rename(columns={"transcript_gene_id": "gene_id"})
    if "fc_label" not in df.columns:
        df["fc_label"] = ""
    return df[["dataset", "gene_id", "transcript_id", "fc_label", "cond1_mult", "cond2_mult"]]


def main():
    args = parse_args()
    out_tx = Path(args.out)
    out_genes = Path(args.out_genes)
    out_tx.parent.mkdir(parents=True, exist_ok=True)
    out_genes.parent.mkdir(parents=True, exist_ok=True)

    frames = [read_design(Path(p)) for p in args.design]
    if not frames:
        pd.DataFrame(columns=TX_FIELDS).to_csv(out_tx, sep="\t", index=False)
        pd.DataFrame(columns=GENE_FIELDS).to_csv(out_genes, sep="\t", index=False)
        return

    df = pd.concat(frames, ignore_index=True)

    df["transcript_direction"] = np.select(
        [df["cond1_mult"] == df["cond2_mult"],
         df["cond1_mult"] > df["cond2_mult"]],
        ["none", "up_in_cond1"],
        default="up_in_cond2",
    )

    # Per-transcript output: only differential transcripts
    (
        df.loc[df["transcript_direction"] != "none", list(TX_FIELDS)]
        .to_csv(out_tx, sep="\t", index=False)
    )

    # Per-gene aggregation: direction counts via get_dummies, gene class via np.select
    direction_dummies = pd.get_dummies(df["transcript_direction"])
    for col in ("none", "up_in_cond1", "up_in_cond2"):
        if col not in direction_dummies.columns:
            direction_dummies[col] = 0

    gene_counts = (
        df[["dataset", "gene_id", "transcript_id"]]
        .join(direction_dummies)
        .groupby(["dataset", "gene_id"], sort=True)
        .agg(
            n_transcripts=("transcript_id", "count"),
            n_up_in_cond1=("up_in_cond1", "sum"),
            n_up_in_cond2=("up_in_cond2", "sum"),
        )
        .assign(n_transcripts_diff=lambda x: x["n_up_in_cond1"] + x["n_up_in_cond2"])
    )

    gene_counts["gene_truth_class"] = np.select(
        [
            gene_counts["n_transcripts_diff"] == 0,
            gene_counts["n_transcripts_diff"] == 1,
            (gene_counts["n_up_in_cond1"] > 0) & (gene_counts["n_up_in_cond2"] > 0),
        ],
        ["none", "single_shift", "isoform_switch"],
        default="co_directional_shift",
    )

    up1_ids = (
        df.loc[df["transcript_direction"] == "up_in_cond1"]
        .groupby(["dataset", "gene_id"])["transcript_id"]
        .agg(";".join)
        .rename("up_in_cond1_transcripts")
    )
    up2_ids = (
        df.loc[df["transcript_direction"] == "up_in_cond2"]
        .groupby(["dataset", "gene_id"])["transcript_id"]
        .agg(";".join)
        .rename("up_in_cond2_transcripts")
    )

    gene_out = (
        gene_counts
        .join(up1_ids, how="left")
        .join(up2_ids, how="left")
        .fillna({"up_in_cond1_transcripts": "", "up_in_cond2_transcripts": ""})
        .reset_index()
    )

    gene_out[list(GENE_FIELDS)].to_csv(out_genes, sep="\t", index=False)


if __name__ == "__main__":
    main()
