# rules/refs.smk
# Canonicalize all reference paths under work/refs/{ds}/...
# Real reference locations from config are symlinked/copied into those paths.
# This avoids using functions in `output:` (not allowed in Snakemake >=8).

from pathlib import Path
import os

# ----- helpers (your original merge kept) -----
def dataset_cfg(ds_name):
    g = config["defaults"].copy()
    ds = config["datasets"][ds_name].get("defaults", {})
    merged = {**g, **ds}
    merged["ref"] = config["datasets"][ds_name]["ref"]
    merged["gtf"] = config["datasets"][ds_name]["gtf"]
    merged["star_index_dir"] = config["datasets"][ds_name].get(
        "star_index_dir",
        str(Path(config["datasets"][ds_name]["ref"]).parent / "STARindex"),
    )
    merged["hisat2_index_prefix"] = config["datasets"][ds_name].get(
        "hisat2_index_prefix",
        str(Path(config["datasets"][ds_name]["ref"]).with_suffix("")),  # e.g., refs/hg38/genome
    )
    return merged

# ===============================
# 1) Canonicalize input references
# ===============================
# We create canonical symlinks (or copies on Windows if needed) so that all downstream rules
# can use static wildcard-based paths for outputs/inputs.

rule link_refs:
    """
    Link (or copy) configured ref & gtf into canonical paths under work/refs/{ds}/
    """
    output:
        ref = "work/refs/{ds}/genome.fa",
        gtf = "work/refs/{ds}/annotation.gtf"
    params:
        src_ref = lambda w: dataset_cfg(w.ds)["ref"],
        src_gtf = lambda w: dataset_cfg(w.ds)["gtf"]
    run:
        import shutil, sys
        out_ref = Path(output.ref)
        out_gtf = Path(output.gtf)
        out_ref.parent.mkdir(parents=True, exist_ok=True)

        def link_or_copy(src, dst):
            src = Path(src)
            dst = Path(dst)
            if dst.exists():
                try:
                    # Refresh symlink/overwrite file
                    if dst.is_symlink() or dst.is_file():
                        dst.unlink()
                except Exception:
                    pass
            try:
                # Prefer symlink where possible
                os.symlink(src, dst)
            except Exception:
                # Fallback to copy (Windows without admin, FS limitations, etc.)
                shutil.copy2(src, dst)

        link_or_copy(params.src_ref, out_ref)
        link_or_copy(params.src_gtf, out_gtf)

# ===============================
# 2) faidx on canonical genome
# ===============================
rule faidx:
    input:
        ref = "work/refs/{ds}/genome.fa"
    output:
        "work/refs/{ds}/genome.fa.fai"
    threads: 1
    conda: "envs/refs.yaml"
    shell:
        r"""
        samtools faidx "{input.ref}"
        """

# ===============================
# 3) STAR genome index (canonical)
# ===============================
rule star_index:
    input:
        ref  = "work/refs/{ds}/genome.fa",
        gtf  = "work/refs/{ds}/annotation.gtf",
        fai  = "work/refs/{ds}/genome.fa.fai"
    output:
        directory("work/refs/{ds}/STARindex")
    threads: 8
    conda: "envs/refs.yaml"
    params:
        # keep your config-driven params, but they don't affect the output path
        sjdbOverhang = lambda w: int(dataset_cfg(w.ds).get("star", {}).get("sjdbOverhang", 100))
    shell:
        r"""
        mkdir -p "{output}"
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir "{output}" \
             --genomeFastaFiles "{input.ref}" \
             --sjdbGTFfile "{input.gtf}" \
             --sjdbOverhang {params.sjdbOverhang}
        """

# ===============================
# 4) HISAT2 index (canonical)
# ===============================
# HISAT2 creates 8 files with a fixed prefix. Use a marker file as the output target.
rule hisat2_index:
    input:
        ref = "work/refs/{ds}/genome.fa",
        fai = "work/refs/{ds}/genome.fa.fai"
    output:
        touch("work/refs/{ds}/hisat2/.done")
    threads: 8
    conda: "envs/refs.yaml"
    params:
        prefix = "work/refs/{ds}/hisat2/genome"
    shell:
        r"""
        mkdir -p "$(dirname "{params.prefix}")"
        hisat2-build -p {threads} "{input.ref}" "{params.prefix}"
        touch "{output}"
        """

# ===============================
# 5) Transcriptome FASTA (canonical)
# ===============================
rule transcripts_fa:
    input:
        ref = "work/refs/{ds}/genome.fa",
        gtf = "work/refs/{ds}/annotation.gtf"
    output:
        "work/refs/{ds}/transcripts.fa"
    conda: "envs/refs.yaml"
    shell:
        r"""
        gffread "{input.gtf}" -g "{input.ref}" -w "{output}"
        """
