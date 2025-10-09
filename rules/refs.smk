# rules/refs.smk
# Builds STAR and HISAT2 indexes per dataset, using paths from config.datasets.<name>
# Requires that genome.fa and annotation.gtf already exist under refs/<genome>/
# (Or point to your own paths in config.yaml.)

from pathlib import Path

# -- helpers (reuse the same merging logic used in your Snakefile) --
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

def star_index_done(ds):
    # STAR creates multiple files; we use genomeParameters.txt as the "done" flag.
    return Path(dataset_cfg(ds)["star_index_dir"]) / "genomeParameters.txt"

def hisat2_index_outputs(ds):
    prefix = dataset_cfg(ds)["hisat2_index_prefix"]
    return [f"{prefix}.{i}.ht2" for i in range(1, 9)]

# --- Samtools faidx (recommended) ---
rule faidx:
    input:
        ref=lambda w: dataset_cfg(w.ds)["ref"]
    output:
        fai=lambda w: dataset_cfg(w.ds)["ref"] + ".fai"
    threads: 1
    conda: "envs/refs.yaml"
    params:
        outdir=lambda w: str(Path(dataset_cfg(w.ds)["ref"]).parent)
    shell:
        r"""
        mkdir -p {params.outdir}
        samtools faidx {input.ref}
        """

# --- STAR genome index ---
rule star_index:
    input:
        ref=lambda w: dataset_cfg(w.ds)["ref"],
        gtf=lambda w: dataset_cfg(w.ds)["gtf"],
        fai=lambda w: dataset_cfg(w.ds)["ref"] + ".fai"
    output:
        done=lambda w: star_index_done(w.ds)
    threads:
        lambda w: int(dataset_cfg(w.ds).get("threads", 8))
    conda: "envs/refs.yaml"
    params:
        outdir=lambda w: dataset_cfg(w.ds)["star_index_dir"],
        sjdbOverhang=lambda w: int(dataset_cfg(w.ds).get("star", {}).get("sjdbOverhang", 100))
    shell:
        r"""
        mkdir -p {params.outdir}
        STAR --runThreadN {threads} \
             --runMode genomeGenerate \
             --genomeDir {params.outdir} \
             --genomeFastaFiles {input.ref} \
             --sjdbGTFfile {input.gtf} \
             --sjdbOverhang {params.sjdbOverhang}
        # mark completion
        test -s {output.done}
        """

# --- HISAT2 index ---
rule hisat2_index:
    input:
        ref=lambda w: dataset_cfg(w.ds)["ref"],
        fai=lambda w: dataset_cfg(w.ds)["ref"] + ".fai"
    output:
        idx=lambda w: hisat2_index_outputs(w.ds)
    threads:
        lambda w: int(dataset_cfg(w.ds).get("threads", 8))
    conda: "envs/refs.yaml"
    params:
        prefix=lambda w: dataset_cfg(w.ds)["hisat2_index_prefix"]
    shell:
        r"""
        mkdir -p $(dirname {params.prefix})
        hisat2-build -p {threads} {input.ref} {params.prefix}
        """

rule transcripts_fa:
    input:
        ref=lambda w: dataset_cfg(w.ds)["ref"],
        gtf=lambda w: dataset_cfg(w.ds)["gtf"]
    output:
        fa=lambda w: str(Path(dataset_cfg(w.ds)["ref"]).parent / "transcripts.fa")
    conda: "envs/refs.yaml"  # ensure gffread is installed there
    shell:
        r"""
        gffread {input.gtf} -g {input.ref} -w {output.fa}
        """
