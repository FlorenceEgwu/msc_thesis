##Purpose

Snakemake-based RNA‑Seq pipeline for the thesis: "Enhanced RNA‑Seq data analysis workflow in molecular biotechnology:
Parameter refinements of mappers for improved accuracy in alternative splicing detection".
Focus on STAR and HISAT2 mapping with rMATS for alternative splicing; 
supports local and Slurm (HPC) execution; uses conda/mamba environments.

##High-level Structure

- Snakefile — main workflow entry point and rules includes.

- config.yaml — central configuration (samples, references, parameters).

-envs/ — conda environment specs (e.g., refs.yaml, star.yaml, hisat2.yaml, rmats.yaml).

- profiles/ — Snakemake profiles: local/ and slurm/.

- rules/ — modular rule files (aligners, QC, post‑processing, splicing).

- scripts/ — helper scripts (pre/post processing, plotting).

- samples.tsv — sample sheet (IDs, groups/conditions, fastq paths or SRR accessions).

- design.tsv — experimental design metadata to complement samples.tsv.

- README.md — quickstart, requirements, and run commands.


##Quickstart

### Requirements
- Mamba/Conda, Snakemake ≥7.32 (recommended), conda-lock (optional)

### Setup
bash
`mamba env create -f envs/refs.yaml`
`mamba env create -f envs/star.yaml`
`mamba env create -f envs/hisat2.yaml`
`mamba env create -f envs/rmats.yaml`

##Local run (example)
`snakemake --profile profiles/local --use-conda -j 4`

##Slurm (HPC)

`snakemake --profile profiles/slurm --use-conda --rerun-incomplete`

##Expected Inputs

- FASTQ files or accession-driven fetch (as defined in samples.tsv).

- Reference genome + annotation (configured via config.yaml and envs/refs.yaml).

- Design matrix (design.tsv) for downstream DE/AS analyses.

##Outputs (typical)

- Aligner-specific BAMs/metrics (STAR/HISAT2).

- Splicing event tables and summaries (rMATS).

- QC reports and intermediate logs.
