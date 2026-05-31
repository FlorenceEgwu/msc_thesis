#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16gb
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=snakemake_dryrun
#SBATCH --output=logs/snakemake/dryrun_%j.log
#SBATCH --error=logs/snakemake/dryrun_%j.err

# ==== SETUP ====
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CORES="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-4}}"
MINIFORGE_ROOT="${HOME}/pl0296-02/project_data/miniforge3"
MINIFORGE_BIN="${MINIFORGE_ROOT}/bin"
JOB_TMP_PARENT="${TMPDIR:-${SLURM_TMPDIR:-${LOCAL_SCRATCH:-/tmp}}}"
JOB_ROOT="${JOB_TMP_PARENT%/}/snakemake-${SLURM_JOB_ID:-dryrun}"

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

# Install rMATS into the same Miniforge if missing.
if [[ ! -x "${MINIFORGE_BIN}/rmats.py" ]]; then
  "${MINIFORGE_BIN}/mamba" install -c conda-forge -c bioconda rmats -y > /dev/null 2>&1
fi

# Install SUPPA2 into the same Miniforge if missing.
if [[ ! -x "${MINIFORGE_BIN}/suppa.py" ]]; then
  "${MINIFORGE_BIN}/mamba" install -c conda-forge -c bioconda suppa -y > /dev/null 2>&1
fi

# ==== DRY-RUN ====
echo "🧪 Starting Snakemake dry run..."
echo "Job ID: ${SLURM_JOB_ID}"
echo "Timestamp: $(date)"
echo "Dry-run cores: ${CORES}"

cd "${REPO_DIR}"

TOTAL_MEM_MB="${SLURM_MEM_PER_NODE:-}"
if [[ -z "${TOTAL_MEM_MB}" && -n "${SLURM_MEM_PER_CPU:-}" && -n "${CORES:-}" ]]; then
  TOTAL_MEM_MB="$((SLURM_MEM_PER_CPU * CORES))"
fi
TOTAL_MEM_MB="${TOTAL_MEM_MB:-16384}"

echo "Memory budget for Snakemake scheduling: ${TOTAL_MEM_MB} MB"
echo ""

"${MINIFORGE_BIN}/snakemake" \
  -s Snakefile \
  --profile profiles/slurm \
  --cores "${CORES}" \
  --resources "mem_mb=${TOTAL_MEM_MB}" \
  --dry-run \
  --rerun-incomplete \
  --keep-going \
  --printshellcmds

echo ""
echo "✅ Dry run completed (no commands were executed)."
echo "Timestamp: $(date)"
