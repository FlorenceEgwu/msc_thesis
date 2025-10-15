# rules/qc.smk — fixed for Snakemake v8+

# If you prefer a single place for OUTDIR, you can keep it as a plain string:
OUTDIR = "results"

rule fastqc:
    input:
        # Build a list of input reads (handles SE or PE gracefully)
        reads=lambda w: (
            [SAMPLE_ROWS[w.sample]["read1"]]
            + ([SAMPLE_ROWS[w.sample]["read2"]] if SAMPLE_ROWS[w.sample].get("read2") else [])
        )
    output:
        "results/fastqc/{sample}.done"
    conda: "envs/qc.yaml"
    threads: 2
    params:
        outdir = f"{OUTDIR}/fastqc"
    shell:
        r"""
        mkdir -p "{params.outdir}"
        fastqc -t {threads} {input.reads} -o "{params.outdir}"
        touch "{output}"
        """

rule multiqc:
    input:
        # Build the list of FastQC done flags for all samples.
        # Using a function here is allowed; it returns a plain list of strings.
        lambda w: expand("results/fastqc/{sample}.done",
                         sample=list(SAMPLE_ROWS.keys()))
    output:
        "results/multiqc/multiqc_report.html"
    conda: "envs/qc.yaml"
    shell:
        r"""
        mkdir -p "{OUTDIR}/multiqc"
        multiqc -o "{OUTDIR}/multiqc" "{OUTDIR}"
        """
