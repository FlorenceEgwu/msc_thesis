# rules/simulate.smk — Snakemake v8+ safe (no functions in output)
# Generates synthetic reads with Polyester per dataset.
# Falls back to FASTA files shipped under reference/ if dataset-specific transcripts aren't set.

import os
from pathlib import Path
from snakemake.exceptions import WorkflowError

# Defaults pulled from config with fallbacks
SIM_CFG = config.get("simulate", {})
SIM_ENABLED = bool(SIM_CFG.get("enabled", False))
SIM_READLEN = int(SIM_CFG.get("readlen", SIM_CFG.get("read_length", 100)))
SIM_FRAGS   = int(SIM_CFG.get("fragments", SIM_CFG.get("reads_per_sample", 100000)))
SIM_SEED    = int(SIM_CFG.get("seed", 12345))

SOFTWARE_DIR = Path(os.path.expanduser("~/pl0296-02/project_data"))

def find_rscript_bin():
    candidates = [
        SOFTWARE_DIR / "miniforge3/bin/Rscript",
        SOFTWARE_DIR / "bin/Rscript",
    ]
    r_root = SOFTWARE_DIR / "R"
    if r_root.exists():
        candidates.extend(sorted(r_root.glob("R-*/bin/Rscript"), reverse=True))
    for path in candidates:
        if path.is_file() and os.access(path, os.X_OK):
            return path
    raise WorkflowError(
        f"Could not find an Rscript binary under {SOFTWARE_DIR}. Update rules/simulate.smk to point to your installation."
    )

RSCRIPT_BIN = find_rscript_bin()

REFERENCE_DIR = Path("reference")
REFERENCE_FASTAS = sorted(
    {p for pattern in ("*.fa", "*.fasta", "*.fa.gz", "*.fna") for p in REFERENCE_DIR.glob(pattern)}
)
DEFAULT_REFERENCE_FASTA = str(next(iter(REFERENCE_FASTAS), ""))

# canonical output dir under repo (can be overridden in params)
def sim_outdir(ds):
    # allow override per-dataset in config["datasets"][ds]["simulate"]["outdir"]
    dsc = config.get("datasets", {}).get(ds, {})
    return dsc.get("simulate", {}).get("outdir", f"data/sim/{ds}")

def sim_dataset_cfg(ds):
    dsc = config.get("datasets", {}).get(ds, {})
    dataset_override = dsc.get("simulate", {})
    global_per_ds = SIM_CFG.get("per_dataset", {}).get(ds, {})
    merged = {}
    merged.update(global_per_ds)
    merged.update(dataset_override)
    return merged

def sim_transcripts(ds):
    cfg = sim_dataset_cfg(ds)
    transcript_path = cfg.get("transcript_fa") or SIM_CFG.get("transcript_fa") or DEFAULT_REFERENCE_FASTA
    if not transcript_path:
        raise WorkflowError(
            f"No transcript FASTA configured for dataset '{ds}' and none found under {REFERENCE_DIR}"
        )
    return transcript_path

def sim_read_length(ds):
    cfg = sim_dataset_cfg(ds)
    return int(cfg.get("read_length", cfg.get("readlen", SIM_READLEN)))

def sim_reads_per_sample(ds):
    cfg = sim_dataset_cfg(ds)
    return int(cfg.get("reads_per_sample", cfg.get("fragments", SIM_FRAGS)))

def sim_n_samples(ds):
    cfg = sim_dataset_cfg(ds)
    return int(cfg.get("n_samples", SIM_CFG.get("n_samples", 1)))

def sim_seed(ds):
    cfg = sim_dataset_cfg(ds)
    return int(cfg.get("seed", SIM_SEED))

def sim_paired_flag(ds):
    cfg = sim_dataset_cfg(ds)
    paired = cfg.get("paired_end", cfg.get("paired", True))
    return "--paired" if paired else ""

def sim_log_path(ds):
    return f"logs/simulate_{ds}.log"

# Use a single marker file to denote completion to avoid listing many fastq outputs
rule simulate_polyester:
    input:
        transcripts = lambda w: sim_transcripts(w.ds)
    output:
        "data/sim/{ds}/done.flag"
    threads: 2
    params:
        outdir   = lambda w: sim_outdir(w.ds),
        read_length = lambda w: sim_read_length(w.ds),
        reads_per_sample = lambda w: sim_reads_per_sample(w.ds),
        n_samples = lambda w: sim_n_samples(w.ds),
        seed     = lambda w: sim_seed(w.ds),
        paired_flag = lambda w: sim_paired_flag(w.ds),
        log      = lambda w: sim_log_path(w.ds),
        rscript  = str(RSCRIPT_BIN)
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{params.outdir}" logs

        # Run the Polyester script (updated earlier to accept --seed; if it also supports --outdir, pass it)
        # If your script does not support --outdir yet, it will default to results/polyester; we will move files after.
        "{params.rscript}" scripts/simulate_polyester.R \
          --transcripts "{input.transcripts}" \
          --outdir "{params.outdir}" \
          --n_samples {params.n_samples} \
          --reads_per_sample {params.reads_per_sample} \
          --read_length {params.read_length} \
          --seed {params.seed} \
          --log "{params.log}" \
          {params.paired_flag}

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
