from pathlib import Path
import os

# Tool paths - use ~ for home directory
STAR_DIR = os.path.expanduser("~/pl0296-02/project_data/STAR-2.7.11b/source")
HISAT2_DIR = os.path.expanduser("~/pl0296-02/project_data/hisat2-2.2.1")
SAMTOOLS_DIR = os.path.expanduser("~/pl0296-02/project_data/samtools-1.19.2")

# output & temp dirs
OUTDIR = "results"
TMPDIR = "work/tmp"

# Helpers
def _sample_id(obj):
    """Normalize either a wildcard object or plain string to a sample ID."""
    if isinstance(obj, str):
        return obj
    if hasattr(obj, "sample"):
        return obj.sample
    raise TypeError(f"Cannot derive sample id from {obj!r}")


def sample_dataset(obj):
    sid = _sample_id(obj)
    return SAMPLE_ROWS[sid].get("dataset", "dataset")


def sample_mapper(obj):
    """Get the mapper type (STAR or HISAT2) for a sample."""
    sid = _sample_id(obj)
    return SAMPLE_ROWS[sid].get("mapper", "HISAT2").upper()

# STAR mapping rule - only triggered for STAR samples
rule map_star:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2", ""),
        star_index = lambda w: STAR_INDEX_TARGET  # noqa: F821 - defined in Snakefile
    output:
        bam = f"{OUTDIR}/mapping/star/{{run}}/{{sample}}.bam",
        bai = f"{OUTDIR}/mapping/star/{{run}}/{{sample}}.bam.bai"
    threads:
        lambda w: sample_threads(w.sample)  # noqa: F821
    params:
        star_dir = STAR_DIR,
        samtools_dir = SAMTOOLS_DIR,
        output_dir = lambda w: sample_output_dir(w.sample),  # noqa: F821
        read_cmd = lambda w: star_read_command(w.sample),    # noqa: F821
        out_samtype = lambda w: star_outsam_arg(w.sample),   # noqa: F821
        twopass = lambda w: star_twopass_mode(w.sample),     # noqa: F821
        mismatch = lambda w: star_mismatch_arg(w.sample),    # noqa: F821
        genome_load = lambda w: sample_star_cfg(w.sample).get("genomeLoad", "NoSharedMemory"),  # noqa: F821
        strand_field = lambda w: star_strand_field(w.sample),  # noqa: F821
        extra = lambda w: sample_star_cfg(w.sample).get("extra", "")  # noqa: F821
    log:
        f"logs/map_star_{{run}}_{{sample}}.log"
    shell:
        r"""
        mkdir -p {TMPDIR} {params.output_dir}
        
        # Prepare read input command
        READ_CMD="{params.read_cmd}"
        READ_IN="--readFilesIn {input.r1}"
        if [ -n "{input.r2}" ]; then
            READ_IN="--readFilesIn {input.r1} {input.r2}"
        fi
        
        {params.star_dir}/STAR \
          --runThreadN {threads} \
          --genomeDir {input.star_index} \
          $READ_IN \
          --readFilesCommand $READ_CMD \
          --outSAMtype {params.out_samtype} \
          --twopassMode {params.twopass} \
          --outFilterMismatchNoverLmax {params.mismatch} \
          --genomeLoad {params.genome_load} \
          --outSAMstrandField {params.strand_field} \
          --outFileNamePrefix {TMPDIR}/{wildcards.sample}.star. \
          {params.extra} \
          > {log} 2>&1
        
        # Move aligned BAM to output
        mv {TMPDIR}/{wildcards.sample}.star.Aligned.sortedByCoord.out.bam {output.bam}
        
        # Index the BAM file
        {params.samtools_dir}/samtools index -@ {threads} {output.bam}
        """

# HISAT2 mapping rule - only triggered for HISAT2 samples
rule map_hisat2:
    input:
        r1 = lambda w: SAMPLE_ROWS[w.sample].get("read1"),
        r2 = lambda w: SAMPLE_ROWS[w.sample].get("read2", ""),
        hisat2_index = lambda w: f"reference/hisat2_index/chr19_4300000-4800000.1.ht2"
    output:
        bam = f"{OUTDIR}/mapping/hisat2/{{run}}/{{sample}}.bam",
        bai = f"{OUTDIR}/mapping/hisat2/{{run}}/{{sample}}.bam.bai"
    threads:
        lambda w: sample_threads(w.sample)  # noqa: F821
    params:
        hisat2_dir = HISAT2_DIR,
        samtools_dir = SAMTOOLS_DIR,
        prefix = "reference/hisat2_index/chr19_4300000-4800000",
        output_dir = lambda w: sample_output_dir(w.sample),  # noqa: F821
        score_min = lambda w: hisat2_score_min_arg(w.sample),  # noqa: F821
        strandness = lambda w: hisat2_strandness(w.sample),    # noqa: F821
        known_splice = lambda w: hisat2_known_splices(w.sample),  # noqa: F821
        extra = lambda w: sample_hisat_cfg(w.sample).get("extra", ""),  # noqa: F821
        dta = lambda w: "1" if sample_hisat_cfg(w.sample).get("dta", True) else "0"  # noqa: F821
    log:
        f"logs/map_hisat2_{{run}}_{{sample}}.log"
    shell:
        r"""
        mkdir -p {params.output_dir}
        
        # Prepare read input command
        READ_IN="-U {input.r1}"
        if [ -n "{input.r2}" ]; then
            READ_IN="-1 {input.r1} -2 {input.r2}"
        fi

        STRAND_OPT=""
        if [ -n "{params.strandness}" ]; then
            STRAND_OPT="--rna-strandness {params.strandness}"
        fi

        KS_OPT=""
        if [ -n "{params.known_splice}" ]; then
            KS_OPT="--known-splicesite-infile {params.known_splice}"
        fi

        DTA_OPT=""
        if [ "{params.dta}" = "1" ]; then
            DTA_OPT="--dta"
        fi
        
        {params.hisat2_dir}/hisat2 \
          -p {threads} \
          $READ_IN \
          -x {params.prefix} \
          --score-min "{params.score_min}" \
          $STRAND_OPT \
          $KS_OPT \
          $DTA_OPT \
          {params.extra} \
          2> {log} \
        | {params.samtools_dir}/samtools sort -@ {threads} -o {output.bam}
        
        # Index the BAM file
        {params.samtools_dir}/samtools index -@ {threads} {output.bam}
        """
