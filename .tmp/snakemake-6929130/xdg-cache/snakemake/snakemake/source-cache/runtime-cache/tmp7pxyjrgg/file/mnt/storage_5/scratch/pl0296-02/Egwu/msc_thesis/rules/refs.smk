# rules/refs.smk

# STAR Genome Index
rule star_index:
    input:
        fasta = REF_FASTA,
        gtf = REF_GTF
    output:
        touch(STAR_INDEX_DONE)
    threads: 8
    resources:
        mem_mb = STAR_INDEX_MEM_MB
    params:
        star_dir = STAR_DIR,
        genome_dir = STAR_INDEX_TARGET,
        sjdbOverhang = 149,
        genomeSAindexNbases = 8
    log:
        "logs/refs/star_index.log"
    shell:
        r"""
        # Check if index already exists and is valid
        if [ -d "{params.genome_dir}" ] && [ -f "{params.genome_dir}/genomeParameters.txt" ]; then
            echo "STAR index already exists at {params.genome_dir}. Skipping indexing..." > {log}
            touch {output}
            exit 0
        fi
        
        mkdir -p {params.genome_dir} $(dirname {output})
        {params.star_dir}/STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {params.genome_dir} \
             --genomeFastaFiles {input.fasta} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang} \
             --genomeSAindexNbases {params.genomeSAindexNbases} \
             > {log} 2>&1
        touch {output}
        """


# HISAT2 Genome Index
rule hisat2_index:
    input:
        fasta = REF_FASTA
    output:
        touch(HISAT2_INDEX_DONE)
    threads: 8
    resources:
        mem_mb = HISAT2_INDEX_MEM_MB
    params:
        hisat2_dir = HISAT2_DIR,
        index_prefix = HISAT2_INDEX_PREFIX
    log:
        "logs/refs/hisat2_index.log"
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
