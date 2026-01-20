#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=snakemake_dryrun
#SBATCH --output=logs/snakemake_dryrun_%j.log
#SBATCH --error=logs/snakemake_dryrun_%j.err

# ==== SETUP ====
set -euo pipefail

# Initialize conda/mamba environment
export PATH="${HOME}/pl0296-02/project_data/miniforge3/bin:$PATH"
export CONDA_PREFIX="${HOME}/pl0296-02/project_data/miniforge3"

# Source conda initialization
eval "$(${HOME}/pl0296-02/project_data/miniforge3/bin/conda shell.bash hook)"
conda activate base

# Ensure snakemake is available
${HOME}/pl0296-02/project_data/miniforge3/bin/mamba install -c conda-forge -c bioconda snakemake -y > /dev/null 2>&1 || true

# Create logs directory
mkdir -p logs

# ==== DRY-RUN ====
echo "🧪 Starting Snakemake dry run..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Timestamp: $(date)"
echo "Dry-run cores: ${SLURM_NTASKS}"
echo ""

cd /mnt/storage_5/scratch/pl0296-02/Egwu/msc_thesis

${HOME}/pl0296-02/project_data/miniforge3/bin/snakemake \
  -s Snakefile \
  --profile profiles/slurm \
  --cores "${SLURM_NTASKS:-4}" \
  --dry-run \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  --verbose

echo ""
echo "✅ Dry run completed (no commands were executed)."
echo "Timestamp: $(date)"
