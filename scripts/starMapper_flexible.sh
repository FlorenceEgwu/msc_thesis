#!/bin/bash

# ==== CONFIGURATION ====
STAR_DIR=${HOME}/pl0296-02/project_data/STAR-2.7.11b/source

# Arguments: dataset_name outFilterMultimapNmax alignIntronMin outdir
DATASET=$1
OUT_FILTER_MULTI=$2
ALIGN_INTRON_MIN=$3
OUT_DIR=$4

if [ -z "$DATASET" ] || [ -z "$OUT_FILTER_MULTI" ] || [ -z "$ALIGN_INTRON_MIN" ] || [ -z "$OUT_DIR" ]; then
    echo "Usage: $0 <dataset_name> <outFilterMultimapNmax> <alignIntronMin> <outdir>"
    exit 1
fi

# Parse config.yaml for dataset
eval $(python3 -c "
import yaml
import sys
with open('config.yaml') as f:
    config = yaml.safe_load(f)
if '$DATASET' not in config['datasets']:
    print('Error: Dataset $DATASET not found in config.yaml')
    sys.exit(1)
dataset = config['datasets']['$DATASET']
print('INPUT_DIR=' + dataset['input_dir'])
print('STAR_INDEX_DIR=' + dataset['star_index_dir'])
print('SAMPLE_PATTERN=' + dataset['sample_name_pattern'])
")

mkdir -p ${OUT_DIR}

# ==== Map reads ====
echo "🚀 Starting mapping process for dataset ${DATASET}..."
for f in ${INPUT_DIR}/sample_*.fasta; do
    SAMPLE=$(basename ${f} .fasta)
    echo "🧩 Processing sample: ${SAMPLE}"
    echo "   Input: ${f}"

    ${STAR_DIR}/STAR \
        --genomeDir ${STAR_INDEX_DIR} \
        --outSAMtype BAM SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif \
        --readFilesIn ${f} \
        --runThreadN 8 \
        --outFilterMultimapNmax ${OUT_FILTER_MULTI} \
        --alignIntronMin ${ALIGN_INTRON_MIN} \
        --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_
done

echo "🎯 Processing complete for dataset ${DATASET}."
# Note: STAR automatically creates sorted BAM files named ${SAMPLE}_Aligned.sortedByCoord.out.bam