# --- Snakefile (fixed header + rule all) ---

import pandas as pd
from pathlib import Path

configfile: "config.yaml"

# Load rule modules FIRST so their symbols/paths exist when we reference them
include: "rules/refs.smk"
include: "rules/mapping.smk"
include: "rules/qc.smk"
include: "rules/splicing.smk"
include: "rules/simulate.smk"

# Samples table
SAMPLES = pd.read_csv("samples.tsv", sep="\t", dtype=str).fillna("")

# Which datasets are used in samples?
DATASETS = sorted(set(SAMPLES["dataset"]))

# Convenience: expand all BAMs using a helper from mapping rules, or define here if needed.
# If your mapping rules already define bam_path(sample_id), keep using that.
# Otherwise, uncomment this local helper to match your mapper outputs.
# def bam_path(sample_id):
#     # e.g., "results/mapping/{sample}.bam" or similar — adjust to your actual layout
#     return f"results/mapping/{sample_id}.bam"

# If your pipeline defines SAMPLE_ROWS dict elsewhere, you can recreate keys here
# assuming a unique 'sample' column:
SAMPLE_IDS = list(SAMPLES["sample"].unique())

# rMATS contrasts, if present in config (e.g., config["as"]["contrasts"] = [["A","B"], ["A","C"]])
CONTRASTS = []
if "as" in config and isinstance(config["as"].get("contrasts"), list):
    CONTRASTS = [tuple(c) for c in config["as"]["contrasts"] if isinstance(c, (list, tuple)) and len(c) == 2]

# Where multiqc report goes (adjust to your true QC output dir)
OUTDIR = Path("results")

def star_targets():
    # STAR index directory per dataset
    return [f"work/refs/{ds}/STARindex" for ds in DATASETS]

def hisat2_targets():
    # HISAT2 index completion marker per dataset (from refs.smk)
    return [f"work/refs/{ds}/hisat2/.done" for ds in DATASETS]

def bam_targets():
    # If you have bam_path(sample_id) defined in mapping rules, use it.
    # Otherwise fall back to a common default (adjust to your actual mapping output).
    try:
        # mapping.smk likely provides bam_path
        return [bam_path(sid) for sid in SAMPLE_IDS]  # noqa: F821 if not defined
    except NameError:
        # Fallback — adjust pattern to your actual mapper outputs
        return [f"results/mapping/{sid}.bam" for sid in SAMPLE_IDS]

def rmats_targets():
    # Only when AS is enabled
    if not config.get("as", {}).get("enabled", False):
        return []
    targets = []
    # Expect a function like rmats_outdir(ds, g1, g2) from splicing.smk.
    # If not present, fall back to a conventional path.
    for ds in DATASETS:
        for (g1, g2) in CONTRASTS:
            try:
                targets.append(str(rmats_outdir(ds, g1, g2) / "MATS_output" / "done.flag"))  # noqa: F821
            except NameError:
                targets.append(f"results/rmats/{ds}/{g1}_vs_{g2}/MATS_output/done.flag")
    return targets

def simulate_targets():
    # Only when simulation is enabled
    if not config.get("simulate", {}).get("enabled", False):
        return []
    t = []
    for ds in DATASETS:
        # If ds_sim_cfg(ds) exists in simulate.smk, prefer it
        try:
            outdir = Path(ds_sim_cfg(ds).get("outdir", f"data/sim/{ds}"))  # noqa: F821
        except NameError:
            outdir = Path(f"data/sim/{ds}")
        t.append(str(outdir / "done.flag"))
    return t

def multiqc_target():
    return [str(OUTDIR / "multiqc" / "multiqc_report.html")] if config.get("qc", {}).get("multiqc", False) else []

rule all:
    input:
        # Reference indexes
        star_targets(),
        hisat2_targets(),
        # Alignments
        bam_targets(),
        # Alternative splicing outputs (when enabled)
        rmats_targets(),
        # Simulation outputs (when enabled)
        simulate_targets(),
        # MultiQC report (when enabled)
        multiqc_target()
