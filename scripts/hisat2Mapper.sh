#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=map_hisat2
#SBATCH --output=logs/map_hisat2_%j.log
#SBATCH --error=logs/map_hisat2_%j.err

# ==== CONFIGURATION ====
HISAT2_DIR=${HOME}/pl0296-02/project_data/hisat2-2.2.1
SAMTOOLS_DIR=${HOME}/pl0296-02/project_data/samtools-1.19.2
data_dir=data/input/fastq_small
out_dir=data/output/bam_hisat2_041125
sufix=".fastq.gz"
FASTA=reference/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
GTF=reference/Sus_scrofa.Sscrofa11.1.104.gtf
INDEX_PREFIX=reference/hisat2_index/Sus_scrofa.Sscrofa11.1

mkdir -p ${out_dir} reference/hisat2_index

# ==== 1. (Optional) Decompress reference files if needed ====
# gunzip -c "${FASTA}.gz" > "${FASTA}"
# gunzip -c "${GTF}.gz" > "${GTF}"

# ==== 2. Build HISAT2 index if not present ====
if [ ! -e "${INDEX_PREFIX}.1.ht2" ]; then
    echo "🧬 Building HISAT2 index..."
    ${HISAT2_DIR}/hisat2-build -p 8 ${FASTA} ${INDEX_PREFIX}
    echo "✅ Index built at ${INDEX_PREFIX}*"
else
    echo "✅ HISAT2 index already exists, skipping build."
fi

# ==== 3. Map reads ====
for f in ${data_dir}/SRR*_2${sufix}; do
{
    SAMPLE=$(basename ${f} _2${sufix})
    PATH_in=${data_dir}/${SAMPLE}

    echo "Processing sample: ${SAMPLE}"
    echo "Reads: ${PATH_in}_1${sufix} and ${PATH_in}_2${sufix}"

    ${HISAT2_DIR}/hisat2 \
        -x ${INDEX_PREFIX} \
        -1 ${PATH_in}_1${sufix} \
        -2 ${PATH_in}_2${sufix} \
        -p 8 \
        --dta \
        --rna-strandness RF \
        | ${SAMTOOLS_DIR}/samtools sort -@ 8 -o ${out_dir}/${SAMPLE}.sorted.bam -

    ${SAMTOOLS_DIR}/samtools index ${out_dir}/${SAMPLE}.sorted.bam
    echo "✅ ${SAMPLE} done."
}
done

echo "🎯 All samples processed with HISAT2."
