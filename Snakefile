import pandas as pd
import yaml
from pathlib import Path

configfile: "config.yaml"

SAMPLES = pd.read_csv("samples.tsv", sep="\t", dtype=str).fillna("")

include: "rules/refs.smk"

# Which datasets are used in samples?
DATASETS = sorted(set(SAMPLES["dataset"]))

# If you want indexes built automatically, add them to the default targets:


def ref_targets():
    wants = []
    for ds in DATASETS:
        wants.append(star_index_done(ds))           # STAR index flag file
        wants += hisat2_index_outputs(ds)          # HISAT2 .ht2 files
    return wants

# Helper: cascade config (global → dataset → per-sample)


def sample_cfg(row):
    g = config["defaults"].copy()
    ds = config["datasets"][row["dataset"]
                            ]["defaults"] if "defaults" in config["datasets"][row["dataset"]] else {}
    merged = {**g, **ds}
    # apply per-sample overrides (flat & dotted keys)
    for col, val in row.items():
        if col in ("sample", "dataset", "read1", "read2", "layout") or val == "":
            continue
        # handle dotted keys like star.twopassMode
        tgt = merged
        parts = col.split(".")
        for p in parts[:-1]:
            if p not in tgt or not isinstance(tgt[p], dict):
                tgt[p] = {}
            tgt = tgt[p]
        tgt[parts[-1]] = val if val != "" else tgt.get(parts[-1], None)
    # inject references
    dsc = config["datasets"][row["dataset"]]
    merged["ref"] = dsc["ref"]
    merged["gtf"] = dsc["gtf"]
    return merged


# Materialize sample table as a dict
SAMPLE_ROWS = {r.sample: r._asdict() for r in SAMPLES.itertuples(index=False)}

# I/O helpers
OUTDIR = Path(config["outdir"])
TMPDIR = Path(config.get("tmpdir", "tmp"))


def bam_path(s):
    return OUTDIR / "bam" / f"{s}.sorted.bam"


rule all:
    input:
        ref_targets(),
        expand(lambda s: bam_path(s), s=SAMPLE_ROWS.keys()),
        ([rmats_outdir(ds, g1, g2) / "MATS_output" / "done.flag" for ds,g1,g2 in CONTRASTS] if config["as"].get("enabled", False) else []),
        ([Path(ds_sim_cfg(ds).get("outdir", f"data/sim/{ds}")) / "sample1_1.fq.gz" for ds in DATASETS] if config.get("simulate",{}).get("enabled", False) else []),
        OUTDIR / "multiqc" / "multiqc_report.html" if config["qc"].get("multiqc", False) else []

include: "rules/mapping.smk"
include: "rules/qc.smk"
include: "rules/splicing.smk"  # enabled when needed
include: "rules/simulate.smk"


