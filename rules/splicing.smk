# rules/splicing.smk — Snakemake v8+ safe rMATS rules (no functions in output)
import pandas as pd
from pathlib import Path

# --- Load design & config ---
DESIGN = pd.read_csv(config["as"]["design_tsv"], sep="\t", dtype=str).fillna("")
READLEN = int(config["as"].get("read_length", 100))
LIBTYPE = config["as"].get("library_type", "unstranded")

# Try to reuse SAMPLES from Snakefile if present; otherwise read samples.tsv
try:
    SAMPLES  # noqa: F821
except NameError:
    SAMPLES = pd.read_csv("samples.tsv", sep="\t", dtype=str).fillna("")

# --- Build contrasts ---
# If you already define config["as"]["contrasts"] = [["A","B"], ...] in Snakefile/config,
# you can use that instead. Otherwise, auto-pick the first 2 conditions per dataset.
CONTRASTS = []
if "as" in config and isinstance(config["as"].get("contrasts"), list):
    for pair in config["as"]["contrasts"]:
        if isinstance(pair, (list, tuple)) and len(pair) == 2:
            for ds in sorted(DESIGN["dataset"].unique()):
                CONTRASTS.append((ds, pair[0], pair[1]))
else:
    for ds, sub in DESIGN.groupby("dataset"):
        conds = list(sub["condition"].unique())
        if len(conds) >= 2:
            CONTRASTS.append((ds, conds[0], conds[1]))

# --- Helpers to collect BAMs produced by mapping rules ---
def bam_path(sample_id: str) -> str:
    # Must match map_star/map_hisat2 outputs; adjust if you changed the layout there.
    return f"results/mapping/{sample_id}.bam"

def samples_for(ds: str, condition: str):
    sub = SAMPLES[(SAMPLES["dataset"] == ds) & (SAMPLES["condition"] == condition)]
    return list(sub["sample"])

def bams_for(ds: str, condition: str):
    return [bam_path(s) for s in samples_for(ds, condition)]

# Expose an expand() helper for Snakefile rule all if you want to target all contrasts
def rmats_done_targets():
    return [
        f"results/rmats/{ds}/{g1}_vs_{g2}/MATS_output/done.flag"
        for (ds, g1, g2) in CONTRASTS
    ]

# --- rMATS rule ---
# NOTE: outputs use static wildcard-based paths (allowed). Dynamic lists are in input/params only.
rule rmats:
    input:
        # canonical refs from refs.smk
        gtf = "work/refs/{ds}/annotation.gtf",
        # BAMs for each group (lists are allowed in input)
        b1_bams = lambda w: bams_for(w.ds, w.g1),
        b2_bams = lambda w: bams_for(w.ds, w.g2)
    output:
        # single marker so Snakemake knows the run completed
        "results/rmats/{ds}/{g1}_vs_{g2}/MATS_output/done.flag"
    threads: lambda w: int(config.get("as", {}).get("threads", 8))
    conda: "envs/rmats.yaml"
    params:
        outdir = "results/rmats/{ds}/{g1}_vs_{g2}",
        tmp    = "work/rmats/{ds}/{g1}_vs_{g2}",
        readlen = READLEN,
        libtype = LIBTYPE
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" "{params.tmp}"

        # rMATS expects text files listing BAMs (one per line).
        # Create them from the input lists.
        B1="{params.tmp}/b1.txt"
        B2="{params.tmp}/b2.txt"
        printf "%s\n" {input.b1_bams} > "$B1"
        printf "%s\n" {input.b2_bams} > "$B2"

        rmats.py \
          --b1 "$B1" \
          --b2 "$B2" \
          --gtf "{input.gtf}" \
          --readLength {params.readlen} \
          --nthread {threads} \
          --libType {params.libtype} \
          --od "{params.outdir}" \
          --tmp "{params.tmp}" \
          --statoff

        # mark completion
        touch "{output}"
        """
