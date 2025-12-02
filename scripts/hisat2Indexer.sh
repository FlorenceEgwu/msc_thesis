#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=8:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --job-name=index_hisat2
#SBATCH --output=logs/index_hisat2_%j.log
#SBATCH --error=logs/index_hisat2_%j.err

# ==== CONFIGURATION ====
HISAT2_DIR=${HOME}/pl0296-02/project_data/hisat2-2.2.1
SAMTOOLS_DIR=${HOME}/pl0296-02/project_data/samtools-1.19.2
FASTA=reference/chr19_4300000-4800000.fa
GTF=reference/chr19_4300000-4800000_fixed.gtf
INDEX_PREFIX=reference/hisat2_index/chr19_4300000-4800000

mkdir -p reference/hisat2_index

# If the FASTA and GTF files are compressed, uncomment the following lines to decompress them
# gunzip -c "${FASTA}.gz" > "${FASTA}"
# gunzip -c "${GTF}.gz" > "${GTF}"

# ==== 2. Build HISAT2 index ====
echo "🧬 Building HISAT2 index..."
${HISAT2_DIR}/hisat2-build -p 8 ${FASTA} ${INDEX_PREFIX}
echo "✅ Index built at ${INDEX_PREFIX}*"

