#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=27
#SBATCH --mem=80gb
#SBATCH --time=168:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=snakemake_pipeline
#SBATCH --output=logs/snakemake_%j.log
#SBATCH --error=logs/snakemake_%j.err

# ==== SETUP ====
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CORES="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-1}}"
MINIFORGE_ROOT="${HOME}/pl0296-02/project_data/miniforge3"
MINIFORGE_BIN="${MINIFORGE_ROOT}/bin"

JOB_TMP_PARENT="${TMPDIR:-${SLURM_TMPDIR:-${LOCAL_SCRATCH:-/tmp}}}"

JOB_ROOT="${JOB_TMP_PARENT%/}/snakemake-${SLURM_JOB_ID:-manual}"

# Keep Python/Snakemake caches off quota-limited home storage.
export TMPDIR="${JOB_ROOT}/tmp"
export XDG_CACHE_HOME="${JOB_ROOT}/xdg-cache"
export PYTHONPYCACHEPREFIX="${JOB_ROOT}/pycache"
export CONDA_PKGS_DIRS="${JOB_ROOT}/conda-pkgs"
mkdir -p "${TMPDIR}" "${XDG_CACHE_HOME}" "${PYTHONPYCACHEPREFIX}" "${CONDA_PKGS_DIRS}" "${REPO_DIR}/logs"

# Initialize conda/mamba environment
export PATH="${MINIFORGE_BIN}:$PATH"
export CONDA_PREFIX="${MINIFORGE_ROOT}"

# Source conda initialization
eval "$(${MINIFORGE_BIN}/conda shell.bash hook)"
conda activate base

# Install snakemake only if it is missing from the existing environment.
if [[ ! -x "${MINIFORGE_BIN}/snakemake" ]]; then
  "${MINIFORGE_BIN}/mamba" install -c conda-forge -c bioconda snakemake -y > /dev/null 2>&1
fi

# ==== RUN SNAKEMAKE PIPELINE ====
echo "🚀 Starting Snakemake RNA-Seq mapping pipeline..."
echo "Job ID: $SLURM_JOB_ID"
echo "Timestamp: $(date)"
echo "Available cores: ${CORES}"
echo ""

cd "${REPO_DIR}"

TOTAL_MEM_MB="${SLURM_MEM_PER_NODE:-}"
if [[ -z "${TOTAL_MEM_MB}" && -n "${SLURM_MEM_PER_CPU:-}" && -n "${CORES:-}" ]]; then
  TOTAL_MEM_MB="$((SLURM_MEM_PER_CPU * CORES))"
fi
TOTAL_MEM_MB="${TOTAL_MEM_MB:-81920}"

echo "Memory budget for Snakemake scheduling: ${TOTAL_MEM_MB} MB"
echo ""

"${MINIFORGE_BIN}/snakemake" \
  -s Snakefile \
  --profile profiles/slurm \
  --cores "${CORES}" \
  --resources "mem_mb=${TOTAL_MEM_MB}" \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds \
  --verbose

  #add -j command to specify number of parallel jobs

echo ""
echo "✅ Pipeline completed successfully!"
echo "Timestamp: $(date)"
