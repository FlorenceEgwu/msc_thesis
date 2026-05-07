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
RMATS_TRUTH_TABLE = RMATS_CFG.get("truth_table", "")
RMATS_MEM_MB = resource_mb("rmats_mem_mb", 32000)
RMATS_SUMMARY_MEM_MB = resource_mb("rmats_summary_mem_mb", 4000)

RMATS_EVENT_FILES = ["SE", "MXE", "A5SS", "A3SS", "RI"]


def rmats_extra_args() -> str:
    args = []
    if RMATS_NOVEL_SS:
        args.append("--novelSS")
    if RMATS_VARIABLE_READ_LENGTH:
        args.append("--variable-read-length")
    if RMATS_PAIRED_STATS:
        args.append("--paired-stats")
    return " ".join(args)


def rmats_case_from_wildcards(wildcards) -> str:
    return f"{wildcards.dataset}/{wildcards.mapper}/{wildcards.param_group}"


rule rmats_bam_lists:
    input:
        group1=lambda wc: rmats_bams(rmats_case_from_wildcards(wc), RMATS_CONDITION_A.get("name", "group1")),
        group2=lambda wc: rmats_bams(rmats_case_from_wildcards(wc), RMATS_CONDITION_B.get("name", "group2")),
    output:
        b1=OUTDIR + "/rmats/{mapper}/{dataset}/{param_group}/b1.txt",
        b2=OUTDIR + "/rmats/{mapper}/{dataset}/{param_group}/b2.txt",
    params:
        group1_name=lambda wc: RMATS_CONDITION_A.get("name", "group1"),
        group2_name=lambda wc: RMATS_CONDITION_B.get("name", "group2"),
    log:
        "logs/rmats/{mapper}/{dataset}/{param_group}/bam_lists.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.b1})" "$(dirname {log})"
        if [[ "{input.group1}" == "" || "{input.group2}" == "" ]]; then
          echo "Missing rMATS replicates for {wildcards.dataset}/{wildcards.mapper}/{wildcards.param_group}" >&2
          echo "{params.group1_name}: {input.group1}" >&2
          echo "{params.group2_name}: {input.group2}" >&2
          exit 1
        fi
        printf '%s\n' "{input.group1}" | tr ' ' ',' > {output.b1}
        printf '%s\n' "{input.group2}" | tr ' ' ',' > {output.b2}
        echo "{params.group1_name}: $(cat {output.b1})" > {log}
        echo "{params.group2_name}: $(cat {output.b2})" >> {log}
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
          --readLength {params.read_length} \
          --libType {params.libtype} \
          --nthread {threads} \
          --tstat {threads} \
          {params.extra} \
          > {log} 2>&1
        touch {output.done}
        """


rule rmats_summarize_case:
    input:
        done=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/rmats.done",
    output:
        significant=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/significant_events.tsv",
        summary=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/summary.tsv",
    params:
        rmats_out=OUTDIR + "/rmats/{dataset}/{mapper}/{param_group}/out",
        truth=RMATS_TRUTH_TABLE,
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
          --truth-table "{params.truth}" \
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
          printf 'dataset\tmapper\tparam_group\tevent_type\ttotal_events\tsignificant_events\tsignificant_genes\ttruth_supported_events\ttruth_supported_genes\n' > {output}
        else
          awk 'FNR==1 && NR!=1 {{next}} {{print}}' {input} > {output}
        fi
        """
