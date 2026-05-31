#!/usr/bin/env python3
"""Emit two SUPPA2-derived AS event tables for one dataset.

--out         per-transcript aggregation: one row per (gene_id, transcript_id)
              with the union of event types it participates in (used by the
              ground-truth GTF table for per-read AS-event stratification).

--out-events  per-event table: one row per SUPPA2 event_id with its
              alternative_transcripts (inclusion form) and total_transcripts.
              This is what summarize_rmats.py joins against to look up
              "which transcript is the inclusion form for an event of this
              type at this gene."
"""

import argparse
import re
from pathlib import Path

import pandas as pd


FIELDS = ["dataset", "gene_id", "transcript_id", "as_event_type", "event_ids", "as_event_source"]
EVENT_FIELDS = ["dataset", "gene_id", "as_event_type", "event_id",
                "alternative_transcripts", "total_transcripts"]
EVENTS = ("SE", "SS", "MXE", "MX", "RI", "FL", "A3SS", "A5SS", "A3", "A5", "AF", "AL")


def event_type(path):
    text = Path(path).name.upper()
    for event in EVENTS:
        if re.search(rf"(^|[_.:;-]){event}($|[_.:;-])", text):
            return {"MX": "MXE", "A3": "A3SS", "A5": "A5SS"}.get(event, event)
    return "unknown"


def read_ioe(path):
    df = pd.read_csv(path, sep="\t")
    total = df["total_transcripts"].fillna("")
    alternative = df.get("alternative_transcripts", pd.Series("", index=df.index)).fillna("")
    transcript_rows = (
        df.assign(
            as_event_type=event_type(path),
            transcript_id=(total + "," + alternative).str.split(r"[,\s]+"),
        )
        .explode("transcript_id")
        .query("transcript_id != ''")
    )
    return transcript_rows[["gene_id", "transcript_id", "as_event_type", "event_id"]]


def read_ioe_events(path):
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["gene_id", "as_event_type", "event_id",
                                     "alternative_transcripts", "total_transcripts"])
    df = df.assign(
        as_event_type=event_type(path),
        alternative_transcripts=df.get("alternative_transcripts",
                                       pd.Series("", index=df.index)).fillna(""),
        total_transcripts=df.get("total_transcripts",
                                 pd.Series("", index=df.index)).fillna(""),
    )
    return df[["gene_id", "as_event_type", "event_id",
               "alternative_transcripts", "total_transcripts"]]


def join_unique(values):
    return ";".join(sorted({str(v) for v in values if pd.notna(v) and str(v)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--ioe", nargs="*", default=[])
    parser.add_argument("--out", required=True,
                        help="Per-transcript aggregated AS event table.")
    parser.add_argument("--out-events", required=True,
                        help="Per-event AS table with inclusion-form transcript IDs.")
    args = parser.parse_args()

    ioe_paths = [Path(ioe_path) for ioe_path in args.ioe if Path(ioe_path).stat().st_size > 0]

    tx_frames = [read_ioe(p) for p in ioe_paths]
    if tx_frames:
        out_tx = (
            pd.concat(tx_frames)
            .groupby(["gene_id", "transcript_id"], as_index=False)
            .agg(as_event_type=("as_event_type", join_unique),
                 event_ids=("event_id", join_unique))
        )
        out_tx.insert(0, "dataset", args.dataset)
        out_tx["as_event_source"] = "SUPPA2"
    else:
        out_tx = pd.DataFrame(columns=FIELDS)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_tx.to_csv(args.out, sep="\t", index=False, columns=FIELDS)

    ev_frames = [read_ioe_events(p) for p in ioe_paths]
    if ev_frames:
        out_ev = pd.concat(ev_frames, ignore_index=True)
        out_ev.insert(0, "dataset", args.dataset)
    else:
        out_ev = pd.DataFrame(columns=EVENT_FIELDS)
    Path(args.out_events).parent.mkdir(parents=True, exist_ok=True)
    out_ev.to_csv(args.out_events, sep="\t", index=False, columns=EVENT_FIELDS)


if __name__ == "__main__":
    main()
