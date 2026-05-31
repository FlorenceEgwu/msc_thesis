# Spec §6: standard mapping QC parsed from each mapper's own log file.
# Per-sample TSV columns are unified across STAR / HISAT2 — see scripts/parse_mapper_logs.py.

def star_log_path(sample: str) -> str:
    return f"{OUTDIR}/mapping/star/{sample_param_group(sample)}/{sample}.star.Log.final.out"


def hisat2_log_path(sample: str) -> str:
    return f"logs/mapping/hisat2/{sample_param_group(sample)}_{sample}.log"


rule mapper_qc_star:
    input:
        bam=lambda wc: bam_path(wc.sample),
    output:
        f"{OUTDIR}/mapper_qc/star/{{run}}/{{sample}}.mapper_qc.tsv"
    threads: 1
    resources:
        mem_mb=1000
    params:
        dataset=lambda wc: sample_cfg(wc.sample).get("dataset", ""),
        log_path=lambda wc: star_log_path(wc.sample),
    log:
        "logs/mapper_qc/star/{run}/{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {PYTHON_BIN} scripts/parse_mapper_logs.py \
          --mapper star \
          --log {params.log_path} \
          --sample {wildcards.sample} \
          --dataset {params.dataset} \
          --run {wildcards.run} \
          --out {output} \
          > {log} 2>&1
        """


rule mapper_qc_hisat2:
    input:
        bam=lambda wc: bam_path(wc.sample),
    output:
        f"{OUTDIR}/mapper_qc/hisat2/{{run}}/{{sample}}.mapper_qc.tsv"
    threads: 1
    resources:
        mem_mb=1000
    params:
        dataset=lambda wc: sample_cfg(wc.sample).get("dataset", ""),
        log_path=lambda wc: hisat2_log_path(wc.sample),
    log:
        "logs/mapper_qc/hisat2/{run}/{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {PYTHON_BIN} scripts/parse_mapper_logs.py \
          --mapper hisat2 \
          --log {params.log_path} \
          --sample {wildcards.sample} \
          --dataset {params.dataset} \
          --run {wildcards.run} \
          --out {output} \
          > {log} 2>&1
        """


rule merge_mapper_qc:
    input:
        MAPPER_QC_TARGETS
    output:
        ALL_MAPPER_QC_TARGET
    threads: 1
    resources:
        mem_mb=4000
    run:
        import pandas as pd
        frames = [pd.read_csv(f, sep="\t") for f in input]
        if frames:
            pd.concat(frames, ignore_index=True).to_csv(output[0], sep="\t", index=False)
        else:
            pd.DataFrame().to_csv(output[0], sep="\t", index=False)
