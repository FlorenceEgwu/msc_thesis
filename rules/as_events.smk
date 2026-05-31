# Alternative-splicing event annotation for the two-transcript simulated datasets.

rule as_selected_gtf:
    input:
        gtf=REF_GTF,
        fasta=lambda wc: as_event_fasta(wc.dataset)
    output:
        OUTDIR + "/as_events/{dataset}/selected.gtf"
    log:
        "logs/as_events/{dataset}/selected_gtf.log"
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        {PYTHON_BIN} scripts/export_selected_gtf.py \
          --gtf {input.gtf} \
          --fasta {input.fasta} \
          --out {output} \
          > {log} 2>&1
        """


rule suppa2_generate_events:
    input:
        gtf=OUTDIR + "/as_events/{dataset}/selected.gtf"
    output:
        suppa_dir=directory(OUTDIR + "/as_events/{dataset}/suppa2")
    params:
        prefix=lambda wc, output: f"{output.suppa_dir}/{wc.dataset}",
        suppa_bin=SUPPA_BIN,
    log:
        "logs/as_events/{dataset}/suppa2_generate_events.log"
    shell:
        r"""
        mkdir -p "{output.suppa_dir}" "$(dirname {log})"
        {params.suppa_bin} generateEvents \
          -i {input.gtf} \
          -o {params.prefix} \
          -f ioe \
          -e SE SS MX RI FL \
          > {log} 2>&1
        """


rule export_as_event_table:
    input:
        suppa_dir=OUTDIR + "/as_events/{dataset}/suppa2"
    output:
        per_transcript=OUTDIR + "/as_events/{dataset}/as_event_table.tsv",
        per_event=OUTDIR + "/as_events/{dataset}/as_event_table_per_event.tsv",
    log:
        "logs/as_events/{dataset}/as_event_table.log"
    shell:
        r"""
        mkdir -p "$(dirname {output.per_transcript})" "$(dirname {log})"
        shopt -s nullglob
        ioe_files=({input.suppa_dir}/*.ioe)
        {PYTHON_BIN} scripts/export_as_event_table.py \
          --dataset {wildcards.dataset} \
          --ioe "${{ioe_files[@]}}" \
          --out {output.per_transcript} \
          --out-events {output.per_event} \
          > {log} 2>&1
        """
