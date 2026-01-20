"""
rules/stringtie.smk

Defines a per-sample StringTie assembly/quant rule. This file expects the
following names to exist in the including Snakefile namespace:

- `bam_for_sample(sample)` -> path to BAM for the given sample
- `STRINGTIE_DIR` -> path to directory containing `stringtie` binary
- `STRINGTIE_GTF` -> annotation GTF path
- `SAMPLE_ROWS` -> dict of sample metadata

The rule writes outputs to `results/stringtie/{sample}.gtf` and
`results/stringtie/{sample}.gene_abund.tab` and logs to `logs/stringtie_{sample}.log`.
"""

import os

STRINGTIE_CONFIG = config.get("stringtie", {}) if "config" in globals() else {}
DEFAULT_STRINGTIE_DIR = os.path.expanduser("~/pl0296-02/project_data/stringtie")

if "STRINGTIE_DIR" not in globals():
    STRINGTIE_DIR = os.path.expanduser(
        STRINGTIE_CONFIG.get("bin_dir")
        or STRINGTIE_CONFIG.get("dir")
        or os.environ.get("STRINGTIE_DIR", DEFAULT_STRINGTIE_DIR)
    )

DEFAULT_STRINGTIE_GTF = "reference/chr19_4300000-4800000_fixed.gtf"


def _stringtie_gtf_from_datasets():
    dataset_cfgs = globals().get("DATASET_CFGS", {})
    sample_rows = globals().get("SAMPLE_ROWS", {})
    ordered_datasets = []

    if sample_rows:
        for row in sample_rows.values():
            dataset = row.get("dataset", "")
            if dataset:
                ordered_datasets.append(dataset)

    if not ordered_datasets and dataset_cfgs:
        ordered_datasets = list(dataset_cfgs.keys())

    for dataset in ordered_datasets:
        gtf = dataset_cfgs.get(dataset, {}).get("gtf")
        if gtf:
            return gtf

    for cfg in dataset_cfgs.values():
        gtf = cfg.get("gtf")
        if gtf:
            return gtf

    return None


if "STRINGTIE_GTF" not in globals():
    STRINGTIE_GTF = (
        STRINGTIE_CONFIG.get("gtf")
        or _stringtie_gtf_from_datasets()
        or DEFAULT_STRINGTIE_GTF
    )

rule stringtie_quant:
    input:
        bam = lambda w: bam_for_sample(w.sample)
    output:
        gtf = "results/stringtie/{sample}.gtf",
        gene = "results/stringtie/{sample}.gene_abund.tab"
    threads:
        lambda w: int(SAMPLE_ROWS[w.sample].get("threads", 8))
    params:
        stringtie_dir = STRINGTIE_DIR,
        gtf = STRINGTIE_GTF
    log:
        "logs/stringtie_{sample}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/stringtie logs

        {params.stringtie_dir}/stringtie {input.bam} \
          -p {threads} -e -G {params.gtf} \
          -A {output.gene} -o {output.gtf} > {log} 2>&1 || true
        """
