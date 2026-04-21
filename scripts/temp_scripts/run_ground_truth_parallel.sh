#!/bin/bash
# Run ground-truth mapping directly for multiple samples in parallel.
#
# Usage:
#   ./run_ground_truth_parallel.sh star_d1_mm10_int21_3 star_d1_mm10_int21_4
#   ./run_ground_truth_parallel.sh --jobs 4 star_d1_mm10_int21_1 star_d1_mm10_int21_2
#
# Supported sample name pattern:
#   <mapper>_d<dataset_number>_<run_bits>_<sample_id>
# Examples:
#   star_d1_mm10_int21_3
#   hisat2_d2_k20_int20_7

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_ground_truth_parallel.sh [--jobs N] <sample1> [sample2 ...]

Options:
  --jobs N     Number of samples to process in parallel (default: 2)
  -h, --help   Show this help message

Examples:
  ./run_ground_truth_parallel.sh star_d1_mm10_int21_3 star_d1_mm10_int21_4
  ./run_ground_truth_parallel.sh --jobs 4 hisat2_d3_k20_int20_1 hisat2_d3_k20_int20_2
EOF
}

JOBS=2
SAMPLES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --jobs)
      [[ $# -ge 2 ]] || { echo "Missing value for --jobs" >&2; exit 1; }
      JOBS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      SAMPLES+=("$1")
      shift
      ;;
  esac
done

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  usage
  exit 1
fi

if ! [[ "${JOBS}" =~ ^[1-9][0-9]*$ ]]; then
  echo "Invalid --jobs value: ${JOBS}" >&2
  exit 1
fi

REPO_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

MINIFORGE_ROOT="${HOME}/pl0296-02/project_data/miniforge3"
MINIFORGE_BIN="${MINIFORGE_ROOT}/bin"
export PATH="${MINIFORGE_BIN}:$PATH"
export CONDA_PREFIX="${MINIFORGE_ROOT}"

eval "$(${MINIFORGE_BIN}/conda shell.bash hook)"
conda activate base

OUTDIR="data/results"
REF_GTF="reference/genome/Mus_musculus.GRCm39.115.gtf"
GTF_TABLE="${OUTDIR}/ground_truth/gtf_exons.tsv"
COMPLEX_EXON_THRESHOLD=2

mkdir -p "${OUTDIR}/ground_truth" "logs/ground_truth"

if [[ ! -f "${GTF_TABLE}" ]]; then
  echo "Ground-truth exon table not found. Building ${GTF_TABLE}..."
  Rscript scripts/export_ground_truth_gtf_table.R \
    --gtf "${REF_GTF}" \
    --complex-exon-threshold "${COMPLEX_EXON_THRESHOLD}" \
    --out "${GTF_TABLE}" \
    > logs/ground_truth/gtf_table.manual.log 2>&1
fi

infer_dataset() {
  local sample="$1"
  if [[ "${sample}" =~ ^[A-Za-z0-9]+_d([0-9]+)_ ]]; then
    printf 'sim_mouse_dataset%s' "${BASH_REMATCH[1]}"
  else
    return 1
  fi
}

infer_truth_reads() {
  local sample="$1"

  if [[ "${sample}" =~ ^[A-Za-z0-9]+_d([0-9]+)_.+_([0-9]+)$ ]]; then
    local dataset_id="${BASH_REMATCH[1]}"
    local sample_id="${BASH_REMATCH[2]}"
    printf 'data/input/sim/polyester_design/dataset%s/sample_%02d_1.fasta\tdata/input/sim/polyester_design/dataset%s/sample_%02d_2.fasta' \
      "${dataset_id}" "${sample_id}" "${dataset_id}" "${sample_id}"
  else
    return 1
  fi
}

run_one_sample() {
  local sample="$1"

  if [[ ! "${sample}" =~ ^([A-Za-z0-9]+)_(.+)_([0-9]+)$ ]]; then
    echo "Skipping invalid sample name '${sample}'. Expected '<mapper>_<run_label>_<sample_id>'." >&2
    return 1
  fi

  local mapper="${BASH_REMATCH[1],,}"
  local run_label="${sample%_*}"
  local dataset
  dataset="$(infer_dataset "${sample}")" || {
    echo "Could not infer dataset from sample '${sample}'." >&2
    return 1
  }

  local truth_paths
  truth_paths="$(infer_truth_reads "${sample}")" || {
    echo "Could not infer truth FASTA paths from sample '${sample}'." >&2
    return 1
  }

  local truth_read1="${truth_paths%%$'\t'*}"
  local truth_read2="${truth_paths##*$'\t'}"

  local bam="${OUTDIR}/mapping/${mapper}/${run_label}/${sample}.bam"
  local outdir="${OUTDIR}/ground_truth/${mapper}/${run_label}"
  local log="logs/ground_truth/${mapper}_${sample}.log"

  if [[ ! -f "${bam}" ]]; then
    echo "Missing BAM for ${sample}: ${bam}" >&2
    return 1
  fi

  if [[ ! -f "${truth_read1}" || ! -f "${truth_read2}" ]]; then
    echo "Missing truth FASTA for ${sample}: ${truth_read1} or ${truth_read2}" >&2
    return 1
  fi

  mkdir -p "${outdir}"

  echo "[$(date '+%F %T')] Starting ${sample}"
  Rscript scripts/mapping_ground_truth.R \
    --bam "${bam}" \
    --gtf-table "${GTF_TABLE}" \
    --truth-read1 "${truth_read1}" \
    --truth-read2 "${truth_read2}" \
    --sample "${sample}" \
    --mapper "${mapper^^}" \
    --dataset "${dataset}" \
    --run "${run_label}" \
    --out-coordinates "${outdir}/${sample}.coordinates.tsv" \
    --out-standard "${outdir}/${sample}.standard_summary.tsv" \
    --out-ground "${outdir}/${sample}.ground_truth_summary.tsv" \
    --out-stratified "${outdir}/${sample}.stratified_summary.tsv" \
    > "${log}" 2>&1
  echo "[$(date '+%F %T')] Finished ${sample}"
}

active_jobs() {
  jobs -rp | wc -l
}

failures=0

for sample in "${SAMPLES[@]}"; do
  while [[ $(active_jobs) -ge ${JOBS} ]]; do
    wait -n || failures=1
  done

  run_one_sample "${sample}" &
done

while [[ $(active_jobs) -gt 0 ]]; do
  wait -n || failures=1
done

if [[ ${failures} -ne 0 ]]; then
  echo "One or more samples failed." >&2
  exit 1
fi

echo "All requested ground-truth jobs completed successfully."
