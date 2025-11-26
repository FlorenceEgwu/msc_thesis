#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=map_star
#SBATCH --output=logs/map_star_%j.log
#SBATCH --error=logs/map_star_%j.err

# ==== CONFIGURATION ====
sciezka_star=${HOME}/pl0296-02/project_data/STAR-2.7.11b/source
data_dir=data/input/fastq_small
out_dir=data/output/bam_star_221125
sufix=".fastq.gz"
INDEX_DIR=reference/star_index_chr19_window

mkdir -p ${out_dir} logs

# ==== Map reads ====
echo "🚀 Starting mapping process..."
for f in ${data_dir}/SRR*_2${sufix}; do
{
    SAMPLE=$(basename ${f} _2${sufix})
    PATH_in=${data_dir}/${SAMPLE}

    echo "🧩 Processing sample: ${SAMPLE}"
    echo "   Reads: ${PATH_in}_1${sufix} and ${PATH_in}_2${sufix}"

    ${sciezka_star}/STAR \
        --genomeDir ${INDEX_DIR} \
        --outSAMtype BAM SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --outSAMstrandField intronMotif \
        --readFilesIn ${PATH_in}_1${sufix} ${PATH_in}_2${sufix} \
        --readFilesCommand zcat \
        --runThreadN 8 \
        --outFileNamePrefix ${out_dir}/${SAMPLE}_
}
done

echo "🎯 Processing complete."
# Note: STAR automatically creates sorted BAM files named ${SAMPLE}_Aligned.sortedByCoord.out.bam