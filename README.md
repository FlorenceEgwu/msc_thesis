## Quickstart

### Requirements
- Mamba/Conda, Snakemake ≥7.32 (recommended), conda-lock (optional)

### Setup
```bash
mamba env create -f envs/refs.yaml
mamba env create -f envs/star.yaml
mamba env create -f envs/hisat2.yaml
mamba env create -f envs/rmats.yaml

##Local run (example)
`snakemake --profile profiles/local --use-conda -j 4`

##Slurm (HPC)

`snakemake --profile profiles/slurm --use-conda --rerun-incomplete`

