#!/bin/bash
# Submit a Snakemake run for one sample only.
#
# Usage:
#   sbatch submit_single_sample.sh hisat2_d3_k20_int20_10
#   sbatch submit_single_sample.sh star_d1_mm10_int21_3
#
# The script infers:
# - mapper: first token before "_"
# - run label: sample name without the trailing "_<sample_number>"

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=snakemake_one_sample_run
#SBATCH --output=logs/snakemake/one_sample_run_%j.log
#SBATCH --error=logs/snakemake/one_sample_run_%j.err

set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: sbatch submit_single_sample.sh <sample_name>"
  echo "Example: sbatch submit_single_sample.sh hisat2_d3_k20_int20_10"
  exit 1
fi

SAMPLE="$1"

if [[ ! "${SAMPLE}" =~ ^([A-Za-z0-9]+)_(.+)_([0-9]+)$ ]]; then
  echo "Sample name '${SAMPLE}' does not match the expected pattern '<mapper>_<run_label>_<sample_id>'."
  exit 1
fi

MAPPER="${BASH_REMATCH[1],,}"
RUN_LABEL="${SAMPLE%_*}"

if [[ "${MAPPER}" != "star" && "${MAPPER}" != "hisat2" ]]; then
  echo "Unsupported mapper '${MAPPER}'. Expected 'star' or 'hisat2'."
  exit 1
fi

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)}"
CORES="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-1}}"
MINIFORGE_ROOT="${HOME}/pl0296-02/project_data/miniforge3"
MINIFORGE_BIN="${MINIFORGE_ROOT}/bin"

JOB_TMP_PARENT="${TMPDIR:-${SLURM_TMPDIR:-${LOCAL_SCRATCH:-/tmp}}}"
JOB_ROOT="${JOB_TMP_PARENT%/}/snakemake-${SLURM_JOB_ID:-manual}"

export TMPDIR="${JOB_ROOT}/tmp"
export XDG_CACHE_HOME="${JOB_ROOT}/xdg-cache"
export PYTHONPYCACHEPREFIX="${JOB_ROOT}/pycache"
export CONDA_PKGS_DIRS="${JOB_ROOT}/conda-pkgs"
mkdir -p "${TMPDIR}" "${XDG_CACHE_HOME}" "${PYTHONPYCACHEPREFIX}" "${CONDA_PKGS_DIRS}" "${REPO_DIR}/logs" "${REPO_DIR}/logs/snakemake"

export PATH="${MINIFORGE_BIN}:$PATH"
export CONDA_PREFIX="${MINIFORGE_ROOT}"

eval "$(${MINIFORGE_BIN}/conda shell.bash hook)"
conda activate base

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

cd "${REPO_DIR}"

OUTDIR="data/results"

TARGETS=(
  "${OUTDIR}/mapping/${MAPPER}/${RUN_LABEL}/${SAMPLE}.bam"
  "${OUTDIR}/mapping/${MAPPER}/${RUN_LABEL}/${SAMPLE}.bam.bai"
  "${OUTDIR}/ground_truth/${MAPPER}/${RUN_LABEL}/${SAMPLE}.coordinates.tsv"
  "${OUTDIR}/ground_truth/${MAPPER}/${RUN_LABEL}/${SAMPLE}.standard_summary.tsv"
  "${OUTDIR}/ground_truth/${MAPPER}/${RUN_LABEL}/${SAMPLE}.ground_truth_summary.tsv"
  "${OUTDIR}/ground_truth/${MAPPER}/${RUN_LABEL}/${SAMPLE}.stratified_summary.tsv"
)

TOTAL_MEM_MB="${SLURM_MEM_PER_NODE:-}"
if [[ -z "${TOTAL_MEM_MB}" && -n "${SLURM_MEM_PER_CPU:-}" && -n "${CORES:-}" ]]; then
  TOTAL_MEM_MB="$((SLURM_MEM_PER_CPU * CORES))"
fi
TOTAL_MEM_MB="${TOTAL_MEM_MB:-32768}"

echo "Starting one-sample Snakemake run"
echo "Job ID: ${SLURM_JOB_ID:-manual}"
echo "Timestamp: $(date)"
echo "Sample: ${SAMPLE}"
echo "Mapper: ${MAPPER}"
echo "Run label: ${RUN_LABEL}"
echo "Available cores: ${CORES}"
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
  "${TARGETS[@]}"

echo ""
echo "One-sample pipeline completed successfully."
echo "Timestamp: $(date)"
