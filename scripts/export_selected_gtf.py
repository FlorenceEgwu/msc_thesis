#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path

import pandas as pd


GTF_COLUMNS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]


def open_text(path):
    return gzip.open(path, "rt") if Path(path).suffix == ".gz" else open(path)


def selected_ids(fasta):
    with open_text(fasta) as handle:
        headers = [line[1:].split()[0].split("|") for line in handle if line.startswith(">")]
    if not headers:
        raise ValueError(f"No transcript IDs found in {fasta}")
    return {h[0] for h in headers}, {h[1] for h in headers if len(h) > 1}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    tx_ids, gene_ids = selected_ids(args.fasta)
    gtf = pd.read_csv(args.gtf, sep="\t", comment="#", names=GTF_COLUMNS, dtype=str)
    tx = gtf["attribute"].str.extract(r'transcript_id "([^"]+)"', expand=False)
    gene = gtf["attribute"].str.extract(r'gene_id "([^"]+)"', expand=False)
    selected = gtf[tx.isin(tx_ids) | ((gtf["feature"] == "gene") & gene.isin(gene_ids))]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    joined = selected[GTF_COLUMNS[0]].astype(str)
    for column in GTF_COLUMNS[1:]:
        joined = joined.str.cat(selected[column].astype(str), sep="\t")
    Path(args.out).write_text(joined.str.cat(sep="\n") + "\n" if len(joined) else "")


if __name__ == "__main__":
    main()
