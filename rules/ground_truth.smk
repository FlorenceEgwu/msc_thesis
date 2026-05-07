rule mapping_ground_truth_per_sample:
    input:
        bam=lambda wc: bam_path(wc.sample),
        bai=lambda wc: bai_path(wc.sample),
        gtf_table=GROUND_TRUTH_GTF_TABLE_TARGET,
        truth_read1=lambda wc: sample_truth_read1(wc.sample),
        truth_read2=lambda wc: sample_truth_read2(wc.sample)
    output:
        coordinates=f"{OUTDIR}/ground_truth/{{mapper}}/{{run}}/{{sample}}.coordinates.tsv",
        standard=f"{OUTDIR}/ground_truth/{{mapper}}/{{run}}/{{sample}}.standard_summary.tsv",
        ground=f"{OUTDIR}/ground_truth/{{mapper}}/{{run}}/{{sample}}.ground_truth_summary.tsv",
        stratified=f"{OUTDIR}/ground_truth/{{mapper}}/{{run}}/{{sample}}.stratified_summary.tsv"
    threads: 1
    resources:
        mem_mb = 8000
    params:
        sample=lambda wc: wc.sample,
        mapper=lambda wc: sample_mapper(wc.sample),
        dataset=lambda wc: sample_cfg(wc.sample).get("dataset", ""),
        run=lambda wc: sample_param_group(wc.sample),
    log:
        f"logs/ground_truth/{{mapper}}/{{run}}/{{sample}}.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.coordinates})" logs/ground_truth

        if [ "{wildcards.mapper}" != "$(printf '%s' "{params.mapper}" | tr '[:upper:]' '[:lower:]')" ]; then
            echo "Mapper mismatch for sample {wildcards.sample}: expected {params.mapper}, got {wildcards.mapper}" > {log}
            exit 1
        fi

        if [ "{wildcards.run}" != "{params.run}" ]; then
            echo "Run label mismatch for sample {wildcards.sample}: expected {params.run}, got {wildcards.run}" > {log}
            exit 1
        fi

        Rscript scripts/mapping_ground_truth.R \
          --bam {input.bam} \
          --gtf-table {input.gtf_table} \
          --truth-read1 {input.truth_read1} \
          --truth-read2 {input.truth_read2} \
          --sample {params.sample} \
          --mapper {params.mapper} \
          --dataset {params.dataset} \
          --run {params.run} \
          --out-coordinates {output.coordinates} \
          --out-standard {output.standard} \
          --out-ground {output.ground} \
          --out-stratified {output.stratified} \
          > {log} 2>&1
        """

rule export_ground_truth_gtf_table:
    input:
        gtf=REF_GTF
    output:
        GROUND_TRUTH_GTF_TABLE_TARGET
    threads: 1
    resources:
        mem_mb = 8000
    params:
        complex_exon_threshold=lambda wc: int(
            config.get("ground_truth", {}).get("complex_exon_threshold", 3)
        )
    log:
        "logs/ground_truth/gtf_table.log"
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        Rscript scripts/export_ground_truth_gtf_table.R \
          --gtf {input.gtf} \
          --complex-exon-threshold {params.complex_exon_threshold} \
          --out {output} \
          > {log} 2>&1
        """

# Merge outputs
#rule merge_ground_truth_coordinates:
#   input:
#        COORDINATE_TARGETS
#    output:
#       ALL_COORDINATES_TARGET
#    threads: 1
#    resources:
#        mem_mb = 8000
#    shell:
#        r"""
#        mkdir -p "$(dirname {output})"
#        awk 'FNR == 1 && NR != 1 {{ next }} {{ print }}' {input} > {output}
#        """

rule merge_ground_truth_standard_summary:
    input:
        STANDARD_SUMMARY_TARGETS
    output:
        ALL_STANDARD_SUMMARY_TARGET
    threads: 1
    resources:
        mem_mb = 8000
    run:
        import pandas as pd
        df = pd.concat([pd.read_csv(f, sep="\t") for f in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule merge_ground_truth_summary:
    input:
        GROUND_TRUTH_SUMMARY_TARGETS
    output:
        ALL_GROUND_TRUTH_SUMMARY_TARGET
    threads: 1
    resources:
        mem_mb = 8000
    run:
        import pandas as pd
        df = pd.concat([pd.read_csv(f, sep="\t") for f in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)

rule merge_ground_truth_stratified:
    input:
        STRATIFIED_SUMMARY_TARGETS
    output:
        ALL_STRATIFIED_SUMMARY_TARGET
    threads: 1
    resources:
        mem_mb = 8000
    run:
        import pandas as pd
        df = pd.concat([pd.read_csv(f, sep="\t") for f in input], ignore_index=True)
        df.to_csv(output[0], sep="\t", index=False)
