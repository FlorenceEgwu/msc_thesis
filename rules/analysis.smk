# Visualization + pivot-table generation from the merged summary TSVs.

rule analysis_report:
    input:
        ground=ALL_GROUND_TRUTH_SUMMARY_TARGET,
        strat=ALL_STRATIFIED_SUMMARY_TARGET,
        standard=ALL_STANDARD_SUMMARY_TARGET,
        mapper_qc=ALL_MAPPER_QC_TARGET,
        rmats=rmats_all_summary_target() if RMATS_ENABLED else [],
    output:
        touch(ANALYSIS_DONE)
    params:
        outdir=f"{OUTDIR}/analysis",
        rmats_arg=lambda wc, input: f"--rmats-summary {input.rmats}" if input.rmats else "",
    resources:
        mem_mb=4000
    log:
        "logs/analysis/analysis_report.log"
    shell:
        r"""
        mkdir -p "{params.outdir}" "$(dirname {log})"
        Rscript scripts/analysis/make_report.R \
          --ground-summary {input.ground} \
          --stratified-summary {input.strat} \
          --standard-summary {input.standard} \
          --mapper-qc-summary {input.mapper_qc} \
          --outdir {params.outdir} \
          {params.rmats_arg} \
          > {log} 2>&1
        """
