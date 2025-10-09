rule fastqc:
    input:
        r1=lambda w: SAMPLE_ROWS[w.sample]["read1"],
        r2=lambda w: SAMPLE_ROWS[w.sample]["read2"]
    output:
        touch(OUTDIR / "fastqc" / "{sample}.done")
    conda: "envs/qc.yaml"
    threads: 2
    shell:
        r"""
        mkdir -p {OUTDIR}/fastqc
        fastqc -t {threads} {input.r1} {input.r2 if input.r2 else ''} -o {OUTDIR}/fastqc
        touch {output}
        """

rule multiqc:
    input:
        expand(OUTDIR / "fastqc" / "{sample}.done", sample=SAMPLE_ROWS.keys())
    output:
        OUTDIR / "multiqc" / "multiqc_report.html"
    conda: "envs/qc.yaml"
    shell:
        r"""
        mkdir -p {OUTDIR}/multiqc
        multiqc -o {OUTDIR}/multiqc {OUTDIR}
        """
