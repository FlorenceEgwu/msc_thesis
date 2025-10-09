# rMATS with two-condition contrasts per dataset.
import pandas as pd
from pathlib import Path

DESIGN = pd.read_csv(config["as"]["design_tsv"], sep="\t", dtype=str).fillna("")
READLEN = int(config["as"].get("read_length", 100))
LIBTYPE = config["as"].get("library_type", "unstranded")

# Build contrasts automatically: for each dataset, take the first two conditions found
CONTRASTS = []
for ds, sub in DESIGN.groupby("dataset"):
    conds = list(sub["condition"].unique())
    if len(conds) >= 2:
        CONTRASTS.append((ds, conds[0], conds[1]))

def bam_for(sample):
    return OUTDIR / "bam" / f"{sample}.sorted.bam"

def rmats_outdir(ds, g1, g2):
    return OUTDIR / "splicing" / "rMATS" / f"{ds}__{g1}_vs_{g2}"

rule rmats_group_lists:
    input:
        bams=lambda w: [bam_for(s) for s in DESIGN.loc[
            (DESIGN.dataset==w.ds) & (DESIGN.condition.isin([w.g1, w.g2]))
        ,"sample"].tolist()]
    output:
        b1=lambda w: rmats_outdir(w.ds, w.g1, w.g2) / "b1.txt",
        b2=lambda w: rmats_outdir(w.ds, w.g1, w.g2) / "b2.txt"
    run:
        dssub = DESIGN[(DESIGN.dataset==wildcards.ds)]
        g1s = dssub[dssub.condition==wildcards.g1]["sample"].tolist()
        g2s = dssub[dssub.condition==wildcards.g2]["sample"].tolist()
        outdir = rmats_outdir(wildcards.ds, wildcards.g1, wildcards.g2)
        outdir.mkdir(parents=True, exist_ok=True)
        (outdir / "b1.txt").write_text(",".join(str(bam_for(s)) for s in g1s))
        (outdir / "b2.txt").write_text(",".join(str(bam_for(s)) for s in g2s))

rule rmats_run:
    input:
        b1=rules.rmats_group_lists.output.b1,
        b2=rules.rmats_group_lists.output.b2,
        gtf=lambda w: config["datasets"][w.ds]["gtf"]
    output:
        touch(rmats_outdir("{ds}", "{g1}", "{g2}") / "MATS_output" / "done.flag")
    threads: 8
    conda: "envs/rmats.yaml"
    params:
        outdir=lambda w: rmats_outdir(w.ds, w.g1, w.g2),
        tmp=lambda w: rmats_outdir(w.ds, w.g1, w.g2) / "tmp",
        readlen=READLEN,
        libtype=LIBTYPE
    shell:
        r"""
        mkdir -p {params.outdir} {params.tmp}
        rmats.py \
          --b1 {input.b1} \
          --b2 {input.b2} \
          --gtf {input.gtf} \
          --readLength {params.readlen} \
          --nthread {threads} \
          --libType {params.libtype} \
          --od {params.outdir} \
          --tmp {params.tmp} \
          --statoff
        touch {output}
        """
