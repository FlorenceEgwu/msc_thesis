#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=24gb
#SBATCH --time=04:00:00
#SBATCH --job-name=stringtie_quant
#SBATCH --output=logs/stringtie_%j.log
#SBATCH --error=logs/stringtie_%j.err

# Usage:
#   ./scripts/stringtieQuant.sh -g reference/annotations.gtf [-i data/output] -o data/output/stringtie -t 8
# Notes:
#   Auto-detects:
#     STAR   BAMs: *_Aligned.sortedByCoord.out.bam
#     HISAT2 BAMs: *.sorted.bam

set -euo pipefail

# ---------- defaults (override with flags) ----------
GTF=reference/chr19_4300000-4800000_fixed.gtf
STAR_BAM_DIR="data/output/bam_star_v2"
HISAT2_BAM_DIR="data/output/bam_hisat2_v2"
OUT_DIR="data/output/stringtie"
THREADS=8
IN_DIR=""   # optional parent/extra dir to scan as well
STRINGTIE_DIR=${HOME}/pl0296-02/project_data/stringtie

usage() {
  cat <<EOF
Usage: $0 -g <annotation.gtf> [-i <scan_dir>] [-o <out_dir>] [-t threads]

Required:
  -g    Reference annotation GTF

Optional:
  -i    Directory to also scan for BAMs (in addition to defaults):
        - STAR dir:   ${STAR_BAM_DIR}
        - HISAT2 dir: ${HISAT2_BAM_DIR}
  -o    Output directory (default: ${OUT_DIR})
  -t    Threads (default: ${THREADS})

Auto-detection patterns:
  STAR:   *_Aligned.sortedByCoord.out.bam
  HISAT2: *.sorted.bam

Outputs per sample:
  - OUT_DIR/<SAMPLE>.gtf
  - OUT_DIR/<SAMPLE>.gene_abund.tab
EOF
  exit 1
}

# --------- parse args ----------
while getopts "g:i:o:t:h" opt; do
  case "$opt" in
    g) GTF="$OPTARG" ;;
    i) IN_DIR="$OPTARG" ;;
    o) OUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${GTF}" ]] && usage
mkdir -p "${OUT_DIR}" logs

echo "GTF:        ${GTF}"
echo "STAR_DIR:   ${STAR_BAM_DIR}"
echo "HISAT2_DIR: ${HISAT2_BAM_DIR}"
[[ -n "${IN_DIR}" ]] && echo "EXTRA DIR:  ${IN_DIR}"
echo "OUT_DIR:    ${OUT_DIR}"
echo "Threads:    ${THREADS}"
echo

# --------- autodetect BAMs ----------
shopt -s nullglob

# Collect candidates from default STAR/HISAT2 dirs
STAR_BAMS=()
HISAT_BAMS=()

# From default STAR dir
for f in "${STAR_BAM_DIR}"/*_Aligned.sortedByCoord.out.bam; do STAR_BAMS+=( "$f" ); done
# From default HISAT2 dir
for f in "${HISAT2_BAM_DIR}"/*.sorted.bam; do HISAT_BAMS+=( "$f" ); done

# If user provided -i, also scan that dir for both patterns
if [[ -n "${IN_DIR}" ]]; then
  for f in "${IN_DIR}"/*_Aligned.sortedByCoord.out.bam; do STAR_BAMS+=( "$f" ); done
  for f in "${IN_DIR}"/*.sorted.bam; do HISAT_BAMS+=( "$f" ); done
fi

# Deduplicate (in case the same file is visible via multiple globs)
declare -A SEEN
DEDUP_STAR=()
for f in "${STAR_BAMS[@]}"; do
  [[ -e "$f" ]] || continue
  if [[ -z "${SEEN[$f]:-}" ]]; then
    SEEN["$f"]=1
    DEDUP_STAR+=( "$f" )
  fi
done
STAR_BAMS=( "${DEDUP_STAR[@]}" )

DEDUP_HISAT=()
for f in "${HISAT_BAMS[@]}"; do
  [[ -e "$f" ]] || continue
  # If a file also matches the STAR pattern (unlikely), keep only once
  if [[ -z "${SEEN[$f]:-}" ]]; then
    SEEN["$f"]=1
    DEDUP_HISAT+=( "$f" )
  fi
done
HISAT_BAMS=( "${DEDUP_HISAT[@]}" )

TOTAL=$(( ${#STAR_BAMS[@]} + ${#HISAT_BAMS[@]} ))
if (( TOTAL == 0 )); then
  echo "No BAMs found. Looked for:"
  echo "  STAR in:   ${STAR_BAM_DIR} ${IN_DIR:+and ${IN_DIR}}"
  echo "  HISAT2 in: ${HISAT2_BAM_DIR} ${IN_DIR:+and ${IN_DIR}}"
  echo "  Patterns:  *_Aligned.sortedByCoord.out.bam  and  *.sorted.bam"
  exit 3
fi

echo "Found ${#STAR_BAMS[@]} STAR BAM(s) and ${#HISAT_BAMS[@]} HISAT2 BAM(s)."
echo

# --------- sample name helpers ----------
sample_from_star() {
  local bam="$1"
  basename "$bam" _Aligned.sortedByCoord.out.bam
}
sample_from_hisat() {
  local bam="$1"
  basename "$bam" .sorted.bam
}

# --------- StringTie runner ----------
run_stringtie() {
  local bam="$1"
  local sample="$2"
  echo "  • ${sample}"
  "${STRINGTIE_DIR}/stringtie" "$bam" \
    -p "${THREADS}" \
    -e -G "${GTF}" \
    -A "${OUT_DIR}/${sample}.gene_abund.tab" \
    -o "${OUT_DIR}/${sample}.gtf"

  # harmless capability probe / placeholder to avoid script failure
  "${STRINGTIE_DIR}/stringtie" --version >/dev/null 2>&1 && \
    "${STRINGTIE_DIR}/stringtie" "$bam" -p "${THREADS}" -e -G "${GTF}" -o /dev/null -A /dev/null \
      -B -C /dev/null >/dev/null 2>&1 || true
}

# --------- process STAR BAMs ----------
if (( ${#STAR_BAMS[@]} > 0 )); then
  echo "Processing STAR BAMs..."
  for BAM in "${STAR_BAMS[@]}"; do
    SAMPLE="$(sample_from_star "$BAM")"
    run_stringtie "$BAM" "$SAMPLE"
  done
  echo
fi

# --------- process HISAT2 BAMs ----------
if (( ${#HISAT_BAMS[@]} > 0 )); then
  echo "Processing HISAT2 BAMs..."
  for BAM in "${HISAT_BAMS[@]}"; do
    SAMPLE="$(sample_from_hisat "$BAM")"
    run_stringtie "$BAM" "$SAMPLE"
  done
  echo
fi

echo "Done. Outputs in: ${OUT_DIR}"
echo "✅ All samples processed with StringTie."
