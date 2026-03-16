# rules/refs.smk
# reference indexing rules based on starIndexer.sh and hisat2Indexer.sh

from pathlib import Path 
import os

# Tool paths - use ~ for home directory
STAR_DIR = os.path.expanduser("~/pl0296-02/project_data/STAR-2.7.11b/source")
HISAT2_DIR = os.path.expanduser("~/pl0296-02/project_data/hisat2-2.2.1")
SAMTOOLS_DIR = os.path.expanduser("~/pl0296-02/project_data/samtools-1.19.2")

# Reference paths
FASTA = "reference/chr19_4300000-4800000.fa"
GTF = "reference/chr19_4300000-4800000_fixed.gtf"
STAR_INDEX_DIR = "reference/star_index_chr19_window"
HISAT2_INDEX_PREFIX = "reference/hisat2_index/chr19_4300000-4800000"

# ===============================
# 1) STAR Genome Index
# ===============================
rule star_index:
    input:
        fasta = FASTA,
        gtf = GTF
    output:
        directory(STAR_INDEX_DIR)
    threads: 8
    params:
        star_dir = STAR_DIR,
        sjdbOverhang = 149,
        genomeSAindexNbases = 8
    log:
        "logs/star_index.log"
    shell:
        r"""
        # Check if index already exists and is valid
        if [ -d "{output}" ] && [ -f "{output}/genomeParameters.txt" ]; then
            echo "STAR index already exists at {output}. Skipping indexing..." > {log}
            exit 0
        fi
        
        mkdir -p {output}
        {params.star_dir}/STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {output} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang} \
             --genomeSAindexNbases {params.genomeSAindexNbases} \
             > {log} 2>&1
        """

# ===============================
# 2) HISAT2 Genome Index
# ===============================
rule hisat2_index:
    input:
        fasta = FASTA
    output:
        touch("logs/hisat2_index.done")
    threads: 8
    params:
        hisat2_dir = HISAT2_DIR,
        index_prefix = HISAT2_INDEX_PREFIX
    log:
        "logs/hisat2_index.log"
    shell:
        r"""
        # Check if index files already exist
        if [ -f "{params.index_prefix}.1.ht2" ] && [ -f "{params.index_prefix}.8.ht2" ]; then
            echo "HISAT2 index already exists at {params.index_prefix}. Skipping indexing..." > {log}
            touch {output}
            exit 0
        fi
        
        mkdir -p $(dirname {params.index_prefix})
        {params.hisat2_dir}/hisat2-build -p {threads} {input.fasta} {params.index_prefix} > {log} 2>&1
        """

