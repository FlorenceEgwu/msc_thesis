## Purpose

Snakemake-based RNA‑Seq pipeline for the thesis: "Enhanced RNA‑Seq data analysis workflow in molecular biotechnology:
Parameter refinements of mappers for improved accuracy in alternative splicing detection".
Focus on STAR and HISAT2 mapping with rMATS for alternative splicing; 
supports local and Slurm (HPC) execution; assumes mapper/QC tools are pre-installed (env specs kept only for reference only).

## High-level Structure

- Snakefile — main workflow entry point and rules includes.

- config.yaml — central configuration (samples, references, parameters).

- profiles/ — Snakemake profiles: local/ and slurm/.

- rules/ — modular rule files (aligners, QC, post‑processing, splicing).

- scripts/ — helper scripts (pre/post processing, plotting).

- README.md — quickstart, requirements, and run commands.


## Quickstart

### Requirements
- Snakemake ≥7.32 plus system-installed STAR, HISAT2, samtools, FastQC, MultiQC, StringTie, and rMATS (conda/mamba optional if you wish to revive the env specs).

### Setup
Ensure the required tools are on your `$PATH` (e.g., via modules or manual installs noted in `rules/refs.smk`).

## Local run (example)
`snakemake --profile profiles/local -j 4`

## Slurm (HPC)

`snakemake --profile profiles/slurm --rerun-incomplete`

## Expected Inputs

- FASTQ files or accession-driven fetch (as defined in samples.tsv).

- Reference genome + annotation (configured via config.yaml and envs/refs.yaml).

- Design matrix (design.tsv) for downstream DE/AS analyses.

## Outputs (typical)

- Aligner-specific BAMs/metrics (STAR/HISAT2).

- Splicing event tables and summaries (rMATS).

- Ground-truth mapping summaries with divisions by read type, transcript complexity,
  mapping-error type, and inferred AS-event class.

- rMATS per-case significant-event tables and merged summaries for simulated
  datasets 2 and 3.

- QC reports and intermediate logs.
