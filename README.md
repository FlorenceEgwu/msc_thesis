# RNA-Seq Snakemake Pipeline

This repository implements a Snakemake-based RNA-Seq workflow for the thesis project "Enhanced RNA-Seq data analysis workflow in molecular biotechnology: Parameter refinements of mappers for improved accuracy in alternative splicing detection." It supports STAR and HISAT2 mapping, ground-truth mapping summaries, and rMATS splicing analysis.

## Repository layout

- `Snakefile` — main workflow entrypoint.
- `config.yaml` — central pipeline configuration.
- `profiles/` — Snakemake profiles for `local` and `slurm` execution.
- `rules/` — modular Snakemake rules for indexing, mapping, ground truth, and rMATS.
- `scripts/` — helper scripts used by rules.
- `data/` — input data, simulated reads, and pipeline outputs.
- `envs/` — reference environment specs for tools.

## Requirements

- Snakemake ≥ 7.32
- STAR
- HISAT2
- samtools
- rMATS
- R with required packages for helper scripts
- Optional: FastQC, MultiQC, StringTie

## Configuration

Primary configuration lives in `config.yaml`.

Important configuration sections:

- `ref`, `gtf` — genome FASTA and annotation GTF.
- `star_index_dir`, `hisat2_index_prefix` — index target locations.
- `outdir`, `tmpdir` — output and temporary directories.
- `datasets` — defines sample templates and parameter groups.
- `tools` — tool directory overrides used by `Snakefile`.
- `rmats` — rMATS pipeline options and sample grouping.

### Notes

### Requirements
- Snakemake ≥7.32 plus system-installed STAR, HISAT2, samtools, FastQC, MultiQC, and StringTie.
- On Slurm, the submit scripts install missing Snakemake and rMATS executables into `~/pl0296-02/project_data/miniforge3`, and `config.yaml` points the rMATS rule at that same Miniforge location.

## Data prerequisites (manual, before any Snakemake run)

The pipeline assumes two manual data-prep stages have already produced their outputs. Neither is wrapped in a Snakemake rule because they depend on Bioconductor packages and on transcript-selection decisions that are outside this repo's scope.

### 1. Selected-transcript FASTAs

Place pre-selected reference transcripts at the following locations. The polyester simulators validate these constraints at startup and abort if anything is off:

- `data/input/selected_transcripts/transcripts_1tx.fa` — exactly **900 transcripts, one per gene**. Order matters: the first 300 must be "easy" (≤2 exons), the next 600 "complex" (>2 exons). FASTA headers must include `transcript_id|gene_id` so the simulator and ground-truth scripts can parse them.
- `data/input/selected_transcripts/transcripts_2tx.fa` — exactly **1800 transcripts in two 900-transcript blocks**. Block 1 is "first-of-pair" (300 easy then 600 complex). Block 2 is "second-of-pair" in the **same gene order** as block 1 (i.e. `gene_id[i] == gene_id[i + 900]`).
- `data/input/selected_transcripts/transcripts_3tx.fa` — same two-block structure as `transcripts_2tx.fa`.

Helper scripts used to build these (run manually):
- `scripts/geneInfoFinder.R` — collects exon counts / transcript metadata from the Ensembl GTF.
- `scripts/geneFastaSelector.R` — picks the transcripts matching the spec rules and writes the FASTA.

The thesis design spec for which transcripts and what ordering is required lives in `resources/settings4sim 4.md`.

### 2. Polyester simulations

Run each simulator manually before invoking Snakemake. They require R with the `polyester`, `Biostrings`, and `dplyr` packages:

```bash
Rscript scripts/polyesterDataset1Simulator.R
Rscript scripts/polyesterDataset2Simulator.R
Rscript scripts/polyesterDataset3Simulator.R
```

Each produces, under `data/input/sim/polyester_design/dataset{N}/`:
- `sample_{NN}_{1,2}.fasta` — 10 paired simulated samples (5 Cond1 + 5 Cond2).
- `transcript_design.tsv` — the canonical per-transcript truth table (used by the rMATS truth rule).
- `sample_design.tsv` — condition / replicate metadata.

These steps live outside the Snakemake DAG because:
- They depend on a fixed seed and only need to be re-run when the simulation design changes.
- Re-running them on every pipeline invocation would be wasteful (heavy IO, no caching benefit).
- The R / Bioconductor environment is independent of the mapping toolchain installed by the submit scripts.

## Running locally

From the repository root:

```bash
snakemake --profile profiles/local -j 4
```

## Running on Slurm

Use the Slurm submission script:

```bash
sbatch submit_pipeline.sh
```

The script will:
- initialize Miniforge in `~/pl0296-02/project_data/miniforge3`
- install `snakemake` and `rMATS` if needed
- run the workflow under the `profiles/slurm` profile

## Inputs and outputs

### Inputs

- Simulated FASTQ files under `data/input/sim/polyester_design/`
- Reference data under `reference/`
- Sample design templates in `config.yaml`

### Outputs

- Aligner-specific BAMs and index files under `data/results/mapping/`
- Ground-truth coordinate and summary tables under `data/results/ground_truth/`
- rMATS outputs and merged summaries under `data/results/rmats/`
- Snakemake logs under `logs/`

## Fixes and validation

The repository has been reviewed for configuration mismatches, and the known `config.yaml` issues were corrected so that STAR parameter groups match their names and HISAT2 group definitions have consistent formatting.
