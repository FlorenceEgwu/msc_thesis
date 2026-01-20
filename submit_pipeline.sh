#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=80gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=snakemake_pipeline
#SBATCH --output=logs/snakemake_%j.log
#SBATCH --error=logs/snakemake_%j.err

# ==== SETUP ====
set -e

# Initialize conda/mamba environment
export PATH="${HOME}/pl0296-02/project_data/miniforge3/bin:$PATH"
export CONDA_PREFIX="${HOME}/pl0296-02/project_data/miniforge3"

# Source conda initialization
eval "$(${HOME}/pl0296-02/project_data/miniforge3/bin/conda shell.bash hook)"
conda activate base

# Install snakemake if not already installed
${HOME}/pl0296-02/project_data/miniforge3/bin/mamba install -c conda-forge -c bioconda snakemake -y > /dev/null 2>&1 || true

# Create logs directory
mkdir -p logs

# ==== RUN SNAKEMAKE PIPELINE ====
echo "🚀 Starting Snakemake RNA-Seq mapping pipeline..."
echo "Job ID: $SLURM_JOB_ID"
echo "Timestamp: $(date)"
echo "Available cores: $SLURM_NTASKS"
echo ""

cd /mnt/storage_3/home/100819/pl0296-02/scratch/Egwu/msc_thesis

${HOME}/pl0296-02/project_data/miniforge3/bin/snakemake \
  -s Snakefile \
  --profile profiles/slurm \
  --cores $SLURM_NTASKS \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  --verbose

echo ""
echo "✅ Pipeline completed successfully!"
echo "Timestamp: $(date)"
