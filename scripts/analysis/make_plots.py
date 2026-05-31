#!/usr/bin/env python3
"""Generate Results-section plots and pivot tables from merged summaries.

Inputs are the four merged TSVs produced by the Snakemake pipeline:
- ground-truth per-sample summary (one row per sample)
- stratified ground-truth summary (read_type x transcript_type x as_event_type x error_type)
- standard mapping summary (BAM-derived QC counts)
- rMATS per-case summary (optional)

Outputs land under --outdir/{plots,tables}/ plus a README.md index.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


AS_DATASETS = ("sim_mouse_dataset2", "sim_mouse_dataset3")
ERROR_TYPES = ("incorrect_position", "incorrect_position_within_region", "incorrect_junction")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ground-summary", required=True)
    p.add_argument("--stratified-summary", required=True)
    p.add_argument("--standard-summary", required=True)
    p.add_argument("--mapper-qc-summary", default="")
    p.add_argument("--rmats-summary", default="")
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def save_fig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def with_label(df):
    return df.assign(case=df["mapper"].str.lower() + ":" + df["run"])


def aggregate_precision(strat, group_cols):
    cols = ["n_reads", "n_correct", "n_incorrect", "n_unmapped"]
    agg = strat.dropna(subset=group_cols).groupby(group_cols, as_index=False)[cols].sum()
    denom = agg["n_correct"] + agg["n_incorrect"]
    agg["precision"] = agg["n_correct"] / denom.where(denom > 0)
    return agg


def grouped_bar(df, *, index, columns, values, title, ylabel, path,
                ylim=None, stacked=False):
    pivot = df.pivot_table(index=index, columns=columns, values=values, aggfunc="mean")
    if pivot.empty:
        return
    pivot = pivot.sort_index()
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(pivot)), 5))
    pivot.plot(kind="bar", ax=ax, stacked=stacked, edgecolor="white", linewidth=0.3)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    if ylim:
        ax.set_ylim(*ylim)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(title=columns, bbox_to_anchor=(1.02, 1), loc="upper left")
    save_fig(fig, path)


def precision_by_structure(strat, plots, tables):
    df = aggregate_precision(strat, ["dataset", "mapper", "run", "transcript_type"])
    df.to_csv(tables / "precision_by_structure.tsv", sep="\t", index=False)
    for dataset, subset in df.groupby("dataset"):
        grouped_bar(with_label(subset), index="case", columns="transcript_type", values="precision",
                    title=f"{dataset} — precision by transcript structure", ylabel="precision",
                    path=plots / f"precision_by_structure_{dataset}.png", ylim=(0, 1))


def precision_by_read_type(strat, plots, tables):
    df = aggregate_precision(strat, ["dataset", "mapper", "run", "read_type"])
    df.to_csv(tables / "precision_by_read_type.tsv", sep="\t", index=False)
    for dataset, subset in df.groupby("dataset"):
        grouped_bar(with_label(subset), index="case", columns="read_type", values="precision",
                    title=f"{dataset} — precision by read type", ylabel="precision",
                    path=plots / f"precision_by_read_type_{dataset}.png", ylim=(0, 1))


def error_breakdown(strat, plots, tables):
    err = strat[strat["error_type"].isin(ERROR_TYPES)]
    if err.empty:
        return
    agg = err.groupby(["dataset", "mapper", "run", "error_type"], as_index=False)["n_incorrect"].sum()
    totals = agg.groupby(["dataset", "mapper", "run"], as_index=False)["n_incorrect"].sum().rename(
        columns={"n_incorrect": "total_incorrect"})
    agg = agg.merge(totals, on=["dataset", "mapper", "run"])
    agg["pct_of_incorrect"] = 100 * agg["n_incorrect"] / agg["total_incorrect"].where(agg["total_incorrect"] > 0)
    agg.to_csv(tables / "error_breakdown.tsv", sep="\t", index=False)
    for dataset, subset in agg.groupby("dataset"):
        grouped_bar(with_label(subset), index="case", columns="error_type", values="pct_of_incorrect",
                    title=f"{dataset} — error-type breakdown",
                    ylabel="% of incorrect alignments",
                    path=plots / f"error_breakdown_{dataset}.png", stacked=True)


def precision_by_as_event(strat, plots, tables):
    eligible = strat[strat["dataset"].isin(AS_DATASETS)]
    if eligible.empty:
        return
    df = aggregate_precision(eligible, ["dataset", "mapper", "run", "as_event_type"])
    df.to_csv(tables / "precision_by_as_event.tsv", sep="\t", index=False)
    for dataset, subset in df.groupby("dataset"):
        grouped_bar(with_label(subset), index="case", columns="as_event_type", values="precision",
                    title=f"{dataset} — precision by AS event", ylabel="precision",
                    path=plots / f"precision_by_as_event_{dataset}.png", ylim=(0, 1))


def standard_qc(standard, plots, tables):
    df = standard.copy()
    df["pct_unmapped"] = 100 * df["unmapped_reads"] / df["total_reads"].where(df["total_reads"] > 0)
    df["pct_multi"] = 100 * df["reads_with_multiple_alignments"] / df["total_reads"].where(df["total_reads"] > 0)
    summary = df.groupby(["dataset", "mapper", "run"], as_index=False)[
        ["pct_mapped", "pct_unmapped", "pct_multi"]
    ].mean()
    summary.to_csv(tables / "standard_qc.tsv", sep="\t", index=False)
    for metric, label in [("pct_mapped", "% mapped"),
                          ("pct_unmapped", "% unmapped"),
                          ("pct_multi", "% multi-aligned")]:
        for dataset, subset in summary.groupby("dataset"):
            grouped_bar(with_label(subset), index="case", columns="mapper", values=metric,
                        title=f"{dataset} — {label}", ylabel=label,
                        path=plots / f"standard_{metric}_{dataset}.png")


def rmats_recovery(rmats, plots, tables):
    if rmats is None or rmats.empty:
        return
    df = rmats.copy()
    df["gene_recovery_rate"] = df["truth_supported_genes"] / df["significant_genes"].where(df["significant_genes"] > 0)
    df["event_recovery_rate"] = df["truth_supported_events"] / df["significant_events"].where(df["significant_events"] > 0)
    df.to_csv(tables / "rmats_truth_recovery.tsv", sep="\t", index=False)
    for dataset, subset in df.groupby("dataset"):
        subset = subset.assign(case=subset["mapper"].str.lower() + ":" + subset["param_group"])
        grouped_bar(subset, index="case", columns="event_type", values="truth_supported_events",
                    title=f"{dataset} — rMATS truth-supported events",
                    ylabel="# truth-supported events",
                    path=plots / f"rmats_truth_events_{dataset}.png")
        grouped_bar(subset, index="case", columns="event_type", values="significant_events",
                    title=f"{dataset} — rMATS significant events",
                    ylabel="# significant events",
                    path=plots / f"rmats_significant_events_{dataset}.png")


def mapper_qc(qc, plots, tables):
    if qc is None or qc.empty:
        return
    metric_cols = [
        ("pct_uniquely_mapped", "% uniquely mapped"),
        ("pct_multi_mapped", "% multi-mapped"),
        ("pct_unmapped", "% unmapped"),
        ("mismatch_rate_pct", "mismatch rate (%)"),
        ("overall_alignment_rate", "HISAT2 overall alignment rate (%)"),
    ]
    available = [m for m in metric_cols if m[0] in qc.columns]
    if not available:
        return
    df = qc.copy()
    for col, _ in available:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    summary = df.groupby(["dataset", "mapper", "run"], as_index=False)[
        [c for c, _ in available]
    ].mean()
    summary.to_csv(tables / "mapper_qc.tsv", sep="\t", index=False)
    for col, label in available:
        for dataset, subset in summary.groupby("dataset"):
            if subset[col].dropna().empty:
                continue
            grouped_bar(with_label(subset), index="case", columns="mapper", values=col,
                        title=f"{dataset} — {label} (from mapper log)",
                        ylabel=label,
                        path=plots / f"mapper_qc_{col}_{dataset}.png")


def rmats_event_class_breakdown(rmats, plots, tables):
    if rmats is None or rmats.empty:
        return
    class_cols = ["events_isoform_switch", "events_single_shift",
                  "events_co_directional", "events_null_gene", "events_unknown_gene"]
    if not all(c in rmats.columns for c in class_cols):
        return
    long = rmats.melt(
        id_vars=["dataset", "mapper", "param_group", "event_type"],
        value_vars=class_cols,
        var_name="event_truth_class", value_name="n_events",
    )
    long.to_csv(tables / "rmats_event_class_breakdown.tsv", sep="\t", index=False)
    per_case = long.groupby(
        ["dataset", "mapper", "param_group", "event_truth_class"], as_index=False
    )["n_events"].sum()
    for dataset, subset in per_case.groupby("dataset"):
        if subset.empty:
            continue
        subset = subset.assign(case=subset["mapper"].str.lower() + ":" + subset["param_group"])
        grouped_bar(subset, index="case", columns="event_truth_class", values="n_events",
                    title=f"{dataset} — rMATS significant events by truth class",
                    ylabel="# significant events",
                    path=plots / f"rmats_event_class_{dataset}.png",
                    stacked=True)


def rmats_suppa2_type_match(rmats, plots, tables):
    if rmats is None or rmats.empty:
        return
    cols = ["events_suppa2_type_present", "events_suppa2_type_absent",
            "events_suppa2_no_annotation"]
    if not all(c in rmats.columns for c in cols):
        return
    long = rmats.melt(
        id_vars=["dataset", "mapper", "param_group", "event_type"],
        value_vars=cols,
        var_name="suppa2_match", value_name="n_events",
    )
    long.to_csv(tables / "rmats_suppa2_type_match.tsv", sep="\t", index=False)
    per_case = long.groupby(
        ["dataset", "mapper", "param_group", "suppa2_match"], as_index=False
    )["n_events"].sum()
    for dataset, subset in per_case.groupby("dataset"):
        if subset.empty:
            continue
        subset = subset.assign(case=subset["mapper"].str.lower() + ":" + subset["param_group"])
        grouped_bar(subset, index="case", columns="suppa2_match", values="n_events",
                    title=f"{dataset} — rMATS events vs SUPPA2 gene-level event types",
                    ylabel="# significant events",
                    path=plots / f"rmats_suppa2_type_match_{dataset}.png",
                    stacked=True)


def rmats_direction_match(rmats, plots, tables):
    if rmats is None or rmats.empty:
        return
    cols = ["events_direction_match", "events_direction_mismatch",
            "events_direction_indeterminate", "events_direction_na"]
    if not all(c in rmats.columns for c in cols):
        return
    long = rmats.melt(
        id_vars=["dataset", "mapper", "param_group", "event_type"],
        value_vars=cols,
        var_name="direction_match", value_name="n_events",
    )
    long.to_csv(tables / "rmats_direction_match.tsv", sep="\t", index=False)
    per_case = long.groupby(
        ["dataset", "mapper", "param_group", "direction_match"], as_index=False
    )["n_events"].sum()
    for dataset, subset in per_case.groupby("dataset"):
        if subset.empty:
            continue
        subset = subset.assign(case=subset["mapper"].str.lower() + ":" + subset["param_group"])
        grouped_bar(subset, index="case", columns="direction_match", values="n_events",
                    title=f"{dataset} — rMATS direction match vs simulation truth",
                    ylabel="# significant events",
                    path=plots / f"rmats_direction_match_{dataset}.png",
                    stacked=True)


def write_index(outdir):
    lines = ["# Pipeline analysis outputs", "",
             "Generated by `scripts/analysis/make_plots.py`.", ""]
    for sub in ("plots", "tables"):
        section_dir = outdir / sub
        if not section_dir.exists():
            continue
        lines.append(f"## {sub.capitalize()}")
        lines.append("")
        for path in sorted(section_dir.glob("*")):
            lines.append(f"- `{sub}/{path.name}`")
        lines.append("")
    (outdir / "README.md").write_text("\n".join(lines))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    plots = outdir / "plots"
    tables = outdir / "tables"
    plots.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    strat = pd.read_csv(args.stratified_summary, sep="\t")
    standard = pd.read_csv(args.standard_summary, sep="\t")
    qc = None
    if args.mapper_qc_summary and Path(args.mapper_qc_summary).exists():
        try:
            qc = pd.read_csv(args.mapper_qc_summary, sep="\t")
        except pd.errors.EmptyDataError:
            qc = None
    rmats = None
    if args.rmats_summary and Path(args.rmats_summary).exists():
        try:
            rmats = pd.read_csv(args.rmats_summary, sep="\t")
        except pd.errors.EmptyDataError:
            rmats = None

    precision_by_structure(strat, plots, tables)
    precision_by_read_type(strat, plots, tables)
    error_breakdown(strat, plots, tables)
    precision_by_as_event(strat, plots, tables)
    standard_qc(standard, plots, tables)
    mapper_qc(qc, plots, tables)
    rmats_recovery(rmats, plots, tables)
    rmats_event_class_breakdown(rmats, plots, tables)
    rmats_suppa2_type_match(rmats, plots, tables)
    rmats_direction_match(rmats, plots, tables)

    write_index(outdir)


if __name__ == "__main__":
    main()
