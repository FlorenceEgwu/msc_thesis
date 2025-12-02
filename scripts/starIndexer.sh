#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=40gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --mail-user=fae.florence1@gmail.com
#SBATCH --account=pl0296-02
#SBATCH --time=06:00:00
#SBATCH --job-name=index_fa
#SBATCH --output=logs/index_fa_%j.log
#SBATCH --error=logs/index_fa_%j.err

sciezka_star=${HOME}/pl0296-02/project_data/STAR-2.7.11b/source
FASTA=reference/chr19_4300000-4800000.fa
GTF=reference/chr19_4300000-4800000_fixed.gtf
INDEX_DIR=reference/star_index_chr19_window
mkdir -p ${INDEX_DIR} logs

# If the FASTA and GTF files are compressed, uncomment the following lines to decompress them
# gunzip -c "${FASTA}.gz" > "${FASTA}"
# gunzip -c "${GTF}.gz" > "${GTF}"

${sciezka_star}/STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir ${INDEX_DIR} \
	--genomeFastaFiles ${FASTA} \
	--sjdbGTFfile ${GTF}\
	--sjdbOverhang 149 \
  	--genomeSAindexNbases 8 // adjust this value based on the genome size



