from pathlib import Path

SIMCFG = config.get("simulate", {})
def ds_sim_cfg(ds):
    return SIMCFG.get("per_dataset", {}).get(ds, {})

def sim_outs(ds):
    c = ds_sim_cfg(ds)
    n = int(c.get("n_samples", 0))
    outdir = Path(c.get("outdir", f"data/sim/{ds}"))
    return [outdir / f"sample{i+1}_1.fasta.gz" for i in range(n)]  # Polyester naming is flexible; we map below

rule build_transcripts_if_needed:
    input:
        lambda w: dataset_cfg(w.ds)["gtf"],
        lambda w: dataset_cfg(w.ds)["ref"]
    output:
        fa=lambda w: Path(dataset_cfg(w.ds)["ref"]).parent / "transcripts.fa"
    conda: "envs/refs.yaml"
    shell:
        r"""
        gffread {input[0]} -g {input[1]} -w {output.fa}
        """

rule simulate_polyester:
    input:
        fa=lambda w: (
            Path(ds_sim_cfg(w.ds).get("transcript_fa","")) if ds_sim_cfg(w.ds).get("transcript_fa")
            else Path(dataset_cfg(w.ds)["ref"]).parent / "transcripts.fa"
        )
    output:
        r1=lambda w: Path(ds_sim_cfg(w.ds).get("outdir", f"data/sim/{w.ds}")) / "sample1_1.fq.gz",
        r2=lambda w: Path(ds_sim_cfg(w.ds).get("outdir", f"data/sim/{w.ds}")) / "sample1_2.fq.gz"
    params:
        outdir=lambda w: ds_sim_cfg(w.ds).get("outdir", f"data/sim/{w.ds}"),
        n=lambda w: int(ds_sim_cfg(w.ds).get("n_samples", 1)),
        reads=lambda w: int(ds_sim_cfg(w.ds).get("reads_per_sample", 1000000)),
        readlen=lambda w: int(ds_sim_cfg(w.ds).get("read_length", 100)),
        paired=lambda w: str(ds_sim_cfg(w.ds).get("paired_end", True)).lower()
    conda: "envs/polyester.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        Rscript scripts/simulate_polyester.R \
          --transcripts {input.fa} \
          --outdir {params.outdir} \
          --n_samples {params.n} \
          --reads_per_sample {params.reads} \
          --read_length {params.readlen} \
          --paired
        # Polyester names files like sampleX_1.fq.gz / sampleX_2.fq.gz by default
        test -s {output.r1}
        """
