# rMATS differential alternative-splicing analysis for simulated datasets 2 and 3.
# One case is dataset x mapper x mapper-parameter group.

RMATS_BIN = os.path.expanduser(TOOL_CFG.get("rmats_bin", "rmats.py"))
PYTHON_BIN = os.path.expanduser(TOOL_CFG.get("python_bin", "python"))
RMATS_THREADS = int(RMATS_CFG.get("threads", 8) or 8)
RMATS_READ_LENGTH = int(RMATS_CFG.get("read_length", 100) or 100)
RMATS_LIB_TYPE = RMATS_CFG.get("lib_type", "fr-unstranded")
RMATS_PAIRED_STATS = str(RMATS_CFG.get("paired_stats", "false")).lower() in {"1", "true", "yes", "y"}
RMATS_FDR = float(RMATS_CFG.get("fdr", 0.05))
RMATS_MIN_ABS_INC_DIFF = float(RMATS_CFG.get("min_abs_inc_level_difference", 0.0))
RMATS_NOVEL_SS = str(RMATS_CFG.get("novel_ss", "false")).lower() in {"1", "true", "yes", "y"}
RMATS_VARIABLE_READ_LENGTH = str(RMATS_CFG.get("variable_read_length", "false")).lower() in {"1", "true", "yes", "y"}
RMATS_TRUTH_TABLE = RMATS_TRUTH_TABLE_TARGET
RMATS_TRUTH_GENES_TABLE = RMATS_CFG.get("truth_genes_table") or f"{OUTDIR}/rmats/simulation_truth_genes.tsv"
RMATS_MEM_MB = resource_mb("rmats_mem_mb", 32000)
RMATS_SUMMARY_MEM_MB = resource_mb("rmats_summary_mem_mb", 4000)


def rmats_extra_args() -> str:
    args = []
    if RMATS_NOVEL_SS:
        args.append("--novelSS")
    if RMATS_VARIABLE_READ_LENGTH:
        args.append("--variable-read-length")
    if RMATS_PAIRED_STATS:
        args.append("--paired-stats")
    return " ".join(args)


if not RMATS_CFG.get("truth_table"):
    rule rmats_simulation_truth:
        input:
            rmats_truth_design_files()
        output:
            transcripts=RMATS_TRUTH_TABLE,
            genes=RMATS_TRUTH_GENES_TABLE
        params:
            design_args=lambda wc, input: " ".join(f"--design {path}" for path in input)
        log:
            "logs/rmats/simulation_truth.log"
        shell:
            r"""
            mkdir -p "$(dirname {output.transcripts})" "$(dirname {log})"
            {PYTHON_BIN} scripts/export_rmats_truth.py \
              {params.design_args} \
              --out {output.transcripts} \
              --out-genes {output.genes} \
              > {log} 2>&1
            """


rule rmats_bam_lists:
    input:
        group1=lambda wc: rmats_bams(wc.dataset, wc.mapper, wc.param_group, RMATS_CONDITION_A),
        group2=lambda wc: rmats_bams(wc.dataset, wc.mapper, wc.param_group, RMATS_CONDITION_B),
    output:
        b1=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/b1.txt",
        b2=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/b2.txt",
    params:
        group1_name=lambda wc: RMATS_CONDITION_A.get("name", "group1"),
        group2_name=lambda wc: RMATS_CONDITION_B.get("name", "group2"),
        group1_bams=lambda wc, input: ",".join(input.group1),
        group2_bams=lambda wc, input: ",".join(input.group2),
    log:
        "logs/rmats/{dataset}/{mapper}/{param_group}/bam_lists.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.b1})" "$(dirname {log})"
        if [[ -z "{params.group1_bams}" || -z "{params.group2_bams}" ]]; then
          echo "Missing rMATS replicates for {wildcards.dataset}/{wildcards.mapper}/{wildcards.param_group}" > {log}
          echo "{params.group1_name}: {params.group1_bams}" >> {log}
          echo "{params.group2_name}: {params.group2_bams}" >> {log}
          exit 1
        fi
        printf '%s\n' "{params.group1_bams}" > {output.b1}
        printf '%s\n' "{params.group2_bams}" > {output.b2}
        echo "{params.group1_name}: {params.group1_bams}" > {log}
        echo "{params.group2_name}: {params.group2_bams}" >> {log}
        """


rule rmats_run:
    input:
        b1=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/b1.txt",
        b2=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/b2.txt",
        gtf=REF_GTF,
    output:
        done=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/rmats.done",
    params:
        od=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/out",
        tmp=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/tmp",
        libtype=RMATS_LIB_TYPE,
        read_length=RMATS_READ_LENGTH,
        extra=lambda wc: rmats_extra_args(),
    threads: RMATS_THREADS
    resources:
        mem_mb=RMATS_MEM_MB
    log:
        "logs/rmats/{dataset}/{mapper}/{param_group}/rmats.log"
    shell:
        r"""
        mkdir -p "{params.od}" "{params.tmp}" "$(dirname {log})"
        {RMATS_BIN} \
          --b1 {input.b1} \
          --b2 {input.b2} \
          --gtf {input.gtf} \
          --od {params.od} \
          --tmp {params.tmp} \
          -t paired \
          --readLength {params.read_length} \
          --libType {params.libtype} \
          --allow-clipping \
          --nthread {threads} \
          --tstat 4 \
          {params.extra} \
          > {log} 2>&1
        touch {output.done}
        """


rule rmats_summarize_case:
    input:
        done=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/rmats.done",
        truth_genes=RMATS_TRUTH_GENES_TABLE,
        truth_tx=RMATS_TRUTH_TABLE,
        suppa2_table=lambda wc: f"{OUTDIR}/as_events/{dataset_short_name(wc.dataset)}/as_event_table.tsv",
        suppa2_events=lambda wc: f"{OUTDIR}/as_events/{dataset_short_name(wc.dataset)}/as_event_table_per_event.tsv",
    output:
        significant=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/significant_events.tsv",
        summary=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/summary.tsv",
    params:
        rmats_out=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/out",
        fdr=RMATS_FDR,
        min_abs_inc_diff=RMATS_MIN_ABS_INC_DIFF,
    resources:
        mem_mb=RMATS_SUMMARY_MEM_MB
    log:
        "logs/rmats/{dataset}/{mapper}/{param_group}/summarize.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        {PYTHON_BIN} scripts/summarize_rmats.py \
          --rmats-out {params.rmats_out} \
          --dataset {wildcards.dataset} \
          --mapper {wildcards.mapper} \
          --param-group {wildcards.param_group} \
          --fdr {params.fdr} \
          --min-abs-inc-diff {params.min_abs_inc_diff} \
          --significant-out {output.significant} \
          --summary-out {output.summary} \
          --truth-table "{input.truth_genes}" \
          --per-transcript-truth "{input.truth_tx}" \
          --suppa2-table "{input.suppa2_table}" \
          --suppa2-events "{input.suppa2_events}" \
          > {log} 2>&1
        """


rule rmats_merge_summaries:
    input:
        rmats_summary_targets()
    output:
        rmats_all_summary_target()
    shell:
        r"""
        mkdir -p "$(dirname {output})"
        if [[ "{input}" == "" ]]; then
          printf 'dataset\tmapper\tparam_group\tevent_type\ttotal_events\tsignificant_events\tsignificant_genes\tevents_isoform_switch\tevents_single_shift\tevents_co_directional\tevents_null_gene\tevents_unknown_gene\ttruth_supported_events\ttruth_supported_genes\tevents_suppa2_type_present\tevents_suppa2_type_absent\tevents_suppa2_no_annotation\tevents_direction_match\tevents_direction_mismatch\tevents_direction_indeterminate\tevents_direction_na\n' > {output}
        else
          awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output}
        fi
        """
