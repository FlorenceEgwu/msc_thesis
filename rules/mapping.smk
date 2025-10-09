from pathlib import Path

localrules:  # convenience for quick local testing; remove on HPC
    map_star, map_hisat2

rule map_reads:
    input:
        r1 = lambda wildcards: SAMPLE_ROWS[wildcards.sample]["read1"],
        r2 = lambda wildcards: SAMPLE_ROWS[wildcards.sample]["read2"],
        ref = lambda w: sample_cfg(SAMPLE_ROWS[w.sample])["ref"],
        gtf = lambda w: sample_cfg(SAMPLE_ROWS[w.sample])["gtf"]
    output:
        bam = lambda w: bam_path(w.sample)
    threads:
        lambda w: int(sample_cfg(SAMPLE_ROWS[w.sample]).get("threads", 8))
    params:
        mapper = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("mapper","STAR").upper(),
        cfg = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]),
        layout = lambda w: SAMPLE_ROWS[w.sample]["layout"] or "PE"
    conda:
        lambda w: f"envs/{'star' if sample_cfg(SAMPLE_ROWS[w.sample]).get('mapper','STAR').upper()=='STAR' else 'hisat2'}.yaml"
    message:
        "Mapping {wildcards.sample} with {params.mapper}"
    shell:
        r"""
        set -euo pipefail
        if [ "{params.mapper}" = "STAR" ]; then
            # placeholder command to record mapper selection; real STAR call is in map_star
            echo "Mapper: STAR"
        fi
        """

rule map_star:
    input:
        r1 = rules.map_reads.input.r1,
        r2 = rules.map_reads.input.r2,
        ref = rules.map_reads.input.ref,
        gtf = rules.map_reads.input.gtf,
        star_idx = lambda w: Path(sample_cfg(SAMPLE_ROWS[w.sample])["star_index_dir"]) / "genomeParameters.txt"
    output:
        bam = rules.map_reads.output.bam
    threads:
        rules.map_reads.threads
    params:
        cfg = rules.map_reads.params.cfg,
        layout = rules.map_reads.params.layout,
        genomeDir = lambda w: sample_cfg(SAMPLE_ROWS[w.sample])["star_index_dir"]
    conda: "envs/star.yaml"
    shell:
                r"""
                set -euo pipefail
                mkdir -p {TMPDIR} {OUTDIR}/bam
                # Ensure input files exist
                test -s {input.r1}
                if [ "{params.layout}" = "PE" ]; then
                        test -s {input.r2}
                        readfiles="{input.r1} {input.r2}"
                else
                        readfiles="{input.r1}"
                fi
                STAR \
                    --runThreadN {threads} \
                    --genomeDir {params.genomeDir} \
                    --readFilesIn $readfiles \
                    --readFilesCommand zcat \
                    --outSAMtype {params.cfg[star][outSAMtype][0]} {params.cfg[star][outSAMtype][1]} \
                    --outFileNamePrefix {TMPDIR}/{wildcards.sample}.
                mv {TMPDIR}/{wildcards.sample}Aligned.sortedByCoord.out.bam {output.bam}
                samtools index {output.bam}
                """

rule map_hisat2:
    input:
        r1 = rules.map_reads.input.r1,
        r2 = rules.map_reads.input.r2,
        ref = rules.map_reads.input.ref,
        gtf = rules.map_reads.input.gtf,
        idx = lambda w: [
            f"{sample_cfg(SAMPLE_ROWS[w.sample])['hisat2_index_prefix']}.{i}.ht2" for i in range(1,9)
        ]
    output:
        bam = rules.map_reads.output.bam
    threads:
        rules.map_reads.threads
    params:
        cfg = rules.map_reads.params.cfg,
        layout = rules.map_reads.params.layout,
        prefix = lambda w: sample_cfg(SAMPLE_ROWS[w.sample])["hisat2_index_prefix"]
    conda: "envs/hisat2.yaml"
    shell:
        r"""
        mkdir -p {TMPDIR} {OUTDIR}/bam
        hisat2 \
          -p {threads} \
          {('-1 ' + input.r1 + ' -2 ' + input.r2) if params.layout=='PE' and input.r2 else ('-U ' + input.r1)} \
          -x {params.prefix} \
          --score-min {params.cfg[hisat2][score_min]} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """