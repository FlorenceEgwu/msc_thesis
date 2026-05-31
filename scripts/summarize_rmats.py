#!/usr/bin/env python3
"""Summarize rMATS output for one (dataset, mapper, param_group) case.
Classifies each significant event along three axes:
  gene_truth:  recovered_isoform_switch | recovered_single_shift |
               false_positive_co_directional | false_positive_null | unknown_gene
  suppa2_type: type_present | type_absent | no_annotation
  direction:   match | mismatch | indeterminate | na
"""

import argparse
import sys
from pathlib import Path
import pandas as pd

EVENT_TYPES = ("SE", "MXE", "A5SS", "A3SS", "RI")
RMATS_SUFFIXES = ("MATS.JC.txt", "MATS.JCEC.txt")
SUMMARY_FIELDS = (
    "dataset", "mapper", "param_group", "event_type", "total_events",
    "significant_events", "significant_genes", "events_isoform_switch", "events_single_shift",
    "events_co_directional", "events_null_gene", "events_unknown_gene",
    "truth_supported_events", "truth_supported_genes",
    "events_suppa2_type_present", "events_suppa2_type_absent", "events_suppa2_no_annotation",
    "events_direction_match", "events_direction_mismatch", "events_direction_indeterminate", "events_direction_na",
)
EVENT_CLASS_BY_GENE_CLASS = {
    "isoform_switch": "recovered_isoform_switch", "single_shift": "recovered_single_shift",
    "co_directional_shift": "false_positive_co_directional", "none": "false_positive_null",
}
TRUTH_SUPPORTED_CLASSES = {"recovered_isoform_switch", "recovered_single_shift"}
SIG_PRELUDE = [
    "dataset", "mapper", "param_group", "event_type", "event_truth_class", "gene_truth_class",
    "suppa2_event_type_present", "suppa2_gene_event_types",
    "suppa2_inclusion_transcripts", "expected_inc_diff_sign", "direction_match",
]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rmats-out", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--mapper", required=True)
    parser.add_argument("--param-group", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-abs-inc-diff", type=float, default=0.0)
    parser.add_argument("--significant-out", required=True)
    parser.add_argument("--summary-out", required=True)
    parser.add_argument("--truth-table", default="", help="Per-gene truth TSV.")
    parser.add_argument("--per-transcript-truth", default="", help="Per-transcript truth TSV.")
    parser.add_argument("--suppa2-table", default="", help="SUPPA2 per-gene AS event table.")
    parser.add_argument("--suppa2-events", default="", help="SUPPA2 per-event table.")
    return parser.parse_args()


def load_tsv(path: str, dataset_filter: str = "") -> pd.DataFrame:
    if not path or not Path(path).exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str)
    if dataset_filter and "dataset" in df.columns:
        df = df[df["dataset"] == dataset_filter]
    return df


def _classify_direction(tx_ids, inc_diff, transcript_truth: dict):
    """Return (expected_sign_str, direction_match_label) for one event."""
    if tx_ids is None:
        return "", "na"
    directions = [transcript_truth.get(tid, "none") for tid in tx_ids]
    has_cond1 = any(d == "up_in_cond1" for d in directions)
    has_cond2 = any(d == "up_in_cond2" for d in directions)
    if has_cond1 and has_cond2:
        return "indeterminate", "indeterminate"
    expected = 1 if has_cond1 else (-1 if has_cond2 else 0)
    if expected == 0:
        return "0", "na"
    if pd.isna(inc_diff) or inc_diff == 0:
        return str(expected), "na"
    return str(expected), "match" if (inc_diff > 0) == (expected > 0) else "mismatch"


def summarize_event(args, event_type, gene_truth, transcript_truth,
                    suppa2_gene_types, suppa2_events):
    rmats_out = Path(args.rmats_out)
    source = next((rmats_out / f"{event_type}.{suffix}" for suffix in RMATS_SUFFIXES
                   if (rmats_out / f"{event_type}.{suffix}").exists()), None)
    case = {"dataset": args.dataset, "mapper": args.mapper,
            "param_group": args.param_group, "event_type": event_type}
    if source is None:
        return pd.DataFrame(), {**dict.fromkeys(SUMMARY_FIELDS, 0), **case}

    df = pd.read_csv(source, sep="\t", dtype=str)
    if not df.empty and ("FDR" not in df.columns or "IncLevelDifference" not in df.columns):
        raise ValueError(f"{source} is missing FDR or IncLevelDifference columns")
    if df.empty:
        return pd.DataFrame(), {**dict.fromkeys(SUMMARY_FIELDS, 0), **case}

    fdr_num = pd.to_numeric(df["FDR"], errors="coerce")
    inc_num = pd.to_numeric(df["IncLevelDifference"], errors="coerce")
    sig = df[(fdr_num <= args.fdr) & (inc_num.abs() >= args.min_abs_inc_diff)].copy()
    sig = sig.assign(**case)
    gene_id_col = sig["GeneID"].fillna("") if "GeneID" in sig.columns else pd.Series("", index=sig.index)
    sig["gene_id"] = gene_id_col
    sig["event_truth_class"] = sig["gene_id"].map(gene_truth).map(EVENT_CLASS_BY_GENE_CLASS).fillna("unknown_gene")
    sig["gene_truth_class"] = sig["gene_id"].map(gene_truth).fillna("")
    gene_type_sets = sig["gene_id"].map(suppa2_gene_types)
    sig["suppa2_event_type_present"] = gene_type_sets.apply(
        lambda types: "no_annotation" if not isinstance(types, set)
        else ("type_present" if event_type in types else "type_absent"))
    sig["suppa2_gene_event_types"] = gene_type_sets.apply(
        lambda types: ";".join(sorted(types)) if isinstance(types, set) else "")
    sig["_inclusion"] = sig["gene_id"].apply(
        lambda gid: (lambda m: m[0] if len(m) == 1 else None)(suppa2_events.get((gid, event_type), [])))
    sig["suppa2_inclusion_transcripts"] = sig["_inclusion"].apply(
        lambda txs: ";".join(txs) if txs is not None else "")
    sig_inc_num = pd.to_numeric(sig["IncLevelDifference"], errors="coerce")
    direction_pairs = [_classify_direction(txs, inc, transcript_truth)
                       for txs, inc in zip(sig["_inclusion"], sig_inc_num)]
    sig["expected_inc_diff_sign"] = [pair[0] for pair in direction_pairs]
    sig["direction_match"] = [pair[1] for pair in direction_pairs]
    truth_counts = sig["event_truth_class"].value_counts()
    suppa2_counts = sig["suppa2_event_type_present"].value_counts()
    dir_counts = sig["direction_match"].value_counts()
    named_genes = sig["gene_id"] != ""
    summary = {
        **case,
        "total_events": len(df), "significant_events": len(sig),
        "significant_genes": sig.loc[named_genes, "gene_id"].nunique(),
        "events_isoform_switch": truth_counts.get("recovered_isoform_switch", 0),
        "events_single_shift": truth_counts.get("recovered_single_shift", 0),
        "events_co_directional": truth_counts.get("false_positive_co_directional", 0),
        "events_null_gene": truth_counts.get("false_positive_null", 0), "events_unknown_gene": truth_counts.get("unknown_gene", 0),
        "truth_supported_events": sig["event_truth_class"].isin(TRUTH_SUPPORTED_CLASSES).sum(),
        "truth_supported_genes": sig.loc[sig["event_truth_class"].isin(TRUTH_SUPPORTED_CLASSES) & named_genes, "gene_id"].nunique(),
        "events_suppa2_type_present": suppa2_counts.get("type_present", 0), "events_suppa2_type_absent": suppa2_counts.get("type_absent", 0),
        "events_suppa2_no_annotation": suppa2_counts.get("no_annotation", 0),
        "events_direction_match": dir_counts.get("match", 0), "events_direction_mismatch": dir_counts.get("mismatch", 0),
        "events_direction_indeterminate": dir_counts.get("indeterminate", 0), "events_direction_na": dir_counts.get("na", 0),
    }
    return sig, summary


def main():
    args = parse_args()
    dataset_short = args.dataset.replace("sim_mouse_", "")

    gene_df = load_tsv(args.truth_table, dataset_short)
    if args.truth_table and gene_df.empty:
        print(
            f"WARNING: truth table is empty after filtering for dataset '{dataset_short}' "
            f"(raw dataset arg: '{args.dataset}'). "
            "All significant events will be classified as unknown_gene. "
            "Check that sim_mouse_ prefix stripping and dataset name match the 'dataset' column in the truth TSV.",
            file=sys.stderr,
        )
    gene_col = next((c for c in ("gene_id", "GeneID") if c in gene_df.columns), None)
    if "gene_truth_class" not in gene_df.columns and not gene_df.empty:
        gene_df["gene_truth_class"] = "single_shift"
    gene_truth = (gene_df.set_index(gene_col)["gene_truth_class"].fillna("single_shift").to_dict()
                  if gene_col else {})

    tx_df = load_tsv(args.per_transcript_truth, dataset_short)
    transcript_truth = (tx_df.set_index("transcript_id")["transcript_direction"].to_dict()
                        if not tx_df.empty and "transcript_id" in tx_df.columns else {})

    suppa2_df = load_tsv(args.suppa2_table)
    if not suppa2_df.empty and {"gene_id", "as_event_type"}.issubset(suppa2_df.columns):
        exploded = suppa2_df.assign(as_event_type=suppa2_df["as_event_type"].str.split(";")).explode("as_event_type")
        suppa2_gene_types = (exploded.assign(as_event_type=exploded["as_event_type"].str.strip())
                             .query("as_event_type != ''").groupby("gene_id")["as_event_type"].apply(set).to_dict())
    else:
        suppa2_gene_types = {}

    ev_df = load_tsv(args.suppa2_events, dataset_short)
    if not ev_df.empty and {"gene_id", "as_event_type"}.issubset(ev_df.columns):
        ev_df["_tx_list"] = ev_df.get("alternative_transcripts", pd.Series("", index=ev_df.index)).fillna("").apply(
            lambda alt: [x.strip() for x in alt.split(",") if x.strip()])
        suppa2_events = ev_df.groupby(["gene_id", "as_event_type"])["_tx_list"].apply(list).to_dict()
    else:
        suppa2_events = {}

    sig_frames, summary_rows = [], []
    for event_type in EVENT_TYPES:
        sig_df, summary = summarize_event(
            args, event_type, gene_truth, transcript_truth, suppa2_gene_types, suppa2_events)
        if not sig_df.empty:
            sig_frames.append(sig_df)
        summary_rows.append(summary)

    sig_out = Path(args.significant_out)
    sig_out.parent.mkdir(parents=True, exist_ok=True)
    if sig_frames:
        combined = pd.concat(sig_frames, ignore_index=True)
        private = {col for col in combined.columns if col.startswith("_")}
        out_cols = [col for col in SIG_PRELUDE + [col for col in combined.columns
                    if col not in SIG_PRELUDE and col not in private] if col in combined.columns]
        combined[out_cols].to_csv(sig_out, sep="\t", index=False)
    else:
        pd.DataFrame(columns=SIG_PRELUDE).to_csv(sig_out, sep="\t", index=False)

    summary_out = Path(args.summary_out)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summary_rows)[list(SUMMARY_FIELDS)].to_csv(summary_out, sep="\t", index=False)


if __name__ == "__main__":
    main()
