# rules/simulate.smk — Snakemake v8+ safe (no functions in output)
# Generates synthetic reads with Polyester per dataset.
# Assumes transcripts FASTA is produced by refs.smk at: work/refs/{ds}/transcripts.fa

from pathlib import Path

# Defaults pulled from config with fallbacks
SIM_CFG = config.get("simulate", {})
SIM_ENABLED = bool(SIM_CFG.get("enabled", False))
SIM_READLEN = int(SIM_CFG.get("readlen", SIM_CFG.get("read_length", 100)))
SIM_FRAGS   = int(SIM_CFG.get("fragments", 100000))
SIM_SEED    = int(SIM_CFG.get("seed", 12345))

# canonical output dir under repo (can be overridden in params)
def sim_outdir(ds):
    # allow override per-dataset in config["datasets"][ds]["simulate"]["outdir"]
    dsc = config.get("datasets", {}).get(ds, {})
    return dsc.get("simulate", {}).get("outdir", f"data/sim/{ds}")

# Use a single marker file to denote completion to avoid listing many fastq outputs
rule simulate_polyester:
    input:
        transcripts = "work/refs/{ds}/transcripts.fa"
    output:
        "data/sim/{ds}/done.flag"
    conda:
        "envs/polyester.yaml"
    threads: 2
    params:
        outdir   = lambda w: sim_outdir(w.ds),
        readlen  = SIM_READLEN,
        fragments= SIM_FRAGS,
        seed     = SIM_SEED
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}"

        # Run the Polyester script (updated earlier to accept --seed; if it also supports --outdir, pass it)
        # If your script does not support --outdir yet, it will default to results/polyester; we will move files after.
        Rscript scripts/simulate_polyester.R \
          --fasta "{input.transcripts}" \
          --readlen {params.readlen} \
          --fragments {params.fragments} \
          --seed {params.seed} || true

        # If the script wrote to results/polyester, move outputs into the canonical outdir.
        if [ -d results/polyester ]; then
          find results/polyester -maxdepth 1 -type f -exec mv -f {{}} "{params.outdir}/" \;
          rmdir results/polyester || true
        fi

        # Mark completion
        touch "{output}"
        """

# Optional: convenience target to expose all dataset simulation flags to rule all
# (If you already create a simulate_targets() helper in Snakefile, you can ignore this.)
def all_sim_done_targets(datasets):
    return [f"data/sim/{ds}/done.flag" for ds in datasets]
