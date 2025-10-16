# rules/mapping.smk — Snakemake v8+ safe (no functions in output)
from pathlib import Path

# Canonical output & temp dirs
OUTDIR = "results"
TMPDIR = "work/tmp"

# Expect these objects to exist (defined in Snakefile or elsewhere you include):
# - SAMPLE_ROWS: dict-like, SAMPLE_ROWS[sample_id] -> row dict with keys: read1, read2, dataset, etc.
# - sample_cfg(row): merges defaults + dataset-specific config, returns dict with keys:
#     - ref, gtf, dataset
#     - star_index_dir / hisat2_index_prefix (but we will use canonical refs from refs.smk)
#     - threads, layout, star params, hisat2 params, etc.

# Helpers
def sample_dataset(w):
    return sample_cfg(SAMPLE_ROWS[w.sample]).get("dataset")

def is_pe(w):
    return sample_cfg(SAMPLE_ROWS[w.sample]).get("layout", "PE").upper() == "PE"

# STAR mapping
rule map_star:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2"),
        star_index = lambda w: f"work/refs/{sample_dataset(w)}/STARindex"
    output:
        bam = f"{OUTDIR}/mapping/{{sample}}.bam"
    threads:
        lambda w: int(sample_cfg(SAMPLE_ROWS[w.sample]).get("threads", 8))
    conda: "envs/star.yaml"
    params:
        layout = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("layout", "PE").upper(),
        # allow extra STAR options via config if you like
        extra  = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("star", {}).get("extra", "")
    shell:
        r"""
        mkdir -p {TMPDIR} {OUTDIR}/mapping
        STAR \
          --runThreadN {threads} \
          --genomeDir "{input.star_index}" \
          {("--readFilesIn " + input.r1 + " " + input.r2) if params.layout=="PE" and input.r2 else ("--readFilesIn " + input.r1)} \
          --readFilesCommand "" \
          --outSAMtype BAM SortedByCoordinate \
          --outFileNamePrefix "{TMPDIR}/{{wildcards.sample}}.star." \
          {params.extra}
        mv "{TMPDIR}/{{wildcards.sample}}.star.Aligned.sortedByCoord.out.bam" "{output.bam}"
        samtools index "{output.bam}"
        """

# HISAT2 mapping
rule map_hisat2:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2"),
        # depend on the hisat2 index "done" marker from refs.smk so ordering is correct
        hisat2_done = lambda w: f"work/refs/{sample_dataset(w)}/hisat2/.done"
    output:
        bam = f"{OUTDIR}/mapping/{{sample}}.bam"
    threads:
        lambda w: int(sample_cfg(SAMPLE_ROWS[w.sample]).get("threads", 8))
    conda: "envs/hisat2.yaml"
    params:
        layout    = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("layout", "PE").upper(),
        # canonical hisat2 prefix built by refs.smk
        prefix    = lambda w: f"work/refs/{sample_dataset(w)}/hisat2/genome",
        # pass specific hisat2 params directly to avoid bracket-templating errors
        score_min = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("hisat2", {}).get("score_min", "L,0,-0.6"),
        extra     = lambda w: sample_cfg(SAMPLE_ROWS[w.sample]).get("hisat2", {}).get("extra", "")
    shell:
        r"""
        mkdir -p {TMPDIR} {OUTDIR}/mapping
        hisat2 \
          -p {threads} \
          {("-1 " + input.r1 + " -2 " + input.r2) if params.layout=="PE" and input.r2 else ("-U " + input.r1)} \
          -x "{params.prefix}" \
          --score-min "{params.score_min}" \
          {params.extra} \
        | samtools sort -@ {threads} -o "{output.bam}"
        samtools index "{output.bam}"
        """

# Convenience umbrella rule if you want a single "map" target that selects mapper from config
rule map_reads:
    input:
        # choose STAR or HISAT2 by touching the corresponding BAM
        # (Snakemake will resolve based on which rule can produce the requested BAM)
        # this is primarily for compatibility if rule names are referenced elsewhere
        bam = f"{OUTDIR}/mapping/{{sample}}.bam"
