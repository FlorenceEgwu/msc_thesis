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
