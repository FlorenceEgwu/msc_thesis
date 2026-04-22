# STAR mapping rule - only triggered for STAR samples
rule map_star:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2", ""),
        star_index_done = STAR_INDEX_DONE
    output:
        bam = mapper_bam_pattern("star"),
        bai = mapper_bai_pattern("star")
    threads: 
        lambda w: sample_threads(w.sample)
    resources:
        mem_mb = lambda w: STAR_MAP_MEM_MB
    params:
        star_dir = STAR_DIR,
        genome_dir = STAR_INDEX_TARGET,
        samtools_dir = SAMTOOLS_DIR,
        output_dir = lambda w: sample_output_dir(w.sample),
        expected_run = lambda w: sample_param_group(w.sample),
        mismatch = lambda w: star_mm_arg(w.sample),
        intronMin = lambda w: star_intron_min_arg(w.sample),
        read_files_command = lambda w: star_read_command(w.sample),
    log:
        f"logs/mapping/star/{{run}}_{{sample}}.log"
    shell: 
        r"""
        mkdir -p {params.output_dir} $(dirname {log})

        if [ "{wildcards.run}" != "{params.expected_run}" ]; then
            echo "Run label mismatch for sample {wildcards.sample}: expected {params.expected_run}, got {wildcards.run}" > {log}
            exit 1
        fi

        {params.star_dir}/STAR \
          --runThreadN {threads} \
          --genomeDir {params.genome_dir} \
          --readFilesIn {input.r1} {input.r2} \
          --readFilesCommand {params.read_files_command} \
          --outSAMstrandField intronMotif \
          --genomeLoad NoSharedMemory \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax {params.mismatch} \
          --alignIntronMin {params.intronMin} \
          --outFileNamePrefix "{params.output_dir}/{wildcards.sample}.star." \
          > {log} 2>&1
        
        # Move aligned BAM to output
        mv {params.output_dir}/{wildcards.sample}.star.Aligned.sortedByCoord.out.bam {output.bam}
        
        # Index the BAM file
        {params.samtools_dir}/samtools index -@ {threads} {output.bam}
        """

# HISAT2 mapping rule - only triggered for HISAT2 samples
rule map_hisat2:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2", ""),
        hisat2_index_done = HISAT2_INDEX_DONE
    output:
        bam = mapper_bam_pattern("hisat2"),
        bai = mapper_bai_pattern("hisat2")
    threads:
        lambda w: sample_threads(w.sample)
    resources:
        mem_mb = lambda w: HISAT2_MAP_MEM_MB
    params:
        hisat2_dir = HISAT2_DIR,
        samtools_dir = SAMTOOLS_DIR,
        prefix = HISAT2_INDEX_PREFIX,
        output_dir = lambda w: sample_output_dir(w.sample),
        expected_run = lambda w: sample_param_group(w.sample),
        read_format = lambda w: hisat2_read_format_flag(w.sample),
        k = lambda w: hisat2_k_arg(w.sample),
        min_intronlen = lambda w: hisat2_min_intronlen_arg(w.sample)
    log:
        f"logs/mapping/hisat2/{{run}}_{{sample}}.log"
    shell:
        r"""
        mkdir -p {params.output_dir} $(dirname {log})

        if [ "{wildcards.run}" != "{params.expected_run}" ]; then
            echo "Run label mismatch for sample {wildcards.sample}: expected {params.expected_run}, got {wildcards.run}" > {log}
            exit 1
        fi
        
        # Prepare read input command
        READ_IN="-U {input.r1}"
        if [ -n "{input.r2}" ]; then
            READ_IN="-1 {input.r1} -2 {input.r2}"
        fi

        {params.hisat2_dir}/hisat2 \
          -p {threads} \
          {params.read_format} \
          $READ_IN \
          -x {params.prefix} \
          -k {params.k} \
          --min-intronlen {params.min_intronlen} \
          2> {log} \
        | {params.samtools_dir}/samtools sort -@ {threads} -o {output.bam}
        
        # Index the BAM file
        {params.samtools_dir}/samtools index -@ {threads} {output.bam}
        """
