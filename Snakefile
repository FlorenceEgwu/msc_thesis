import os
import workflowSampleDesign

# Match execution environment
shell.executable("/bin/bash")
shell.prefix("set -euo pipefail")

# ---- Config inputs ----
configfile: "config.yaml"

#ouptut and temp directories
OUTDIR = config.get("outdir", "data/results")
TMPDIR = config.get("tmpdir", "tmp")

# Reference files and index locations
REF_FASTA = config.get("ref", "reference/genome/Mus_musculus.GRCm39.dna.primary_assembly.fa")
REF_GTF = config.get("gtf", "reference/genome/Mus_musculus.GRCm39.115.gtf")
STAR_INDEX_TARGET = config.get("star_index_dir",  "reference/star_index")
STAR_INDEX_DONE = "logs/star_index.done"
HISAT2_INDEX_PREFIX = config.get("hisat2_index_prefix", "reference/hisat2_index/mus_musculus_index")
HISAT2_INDEX_DONE = "logs/hisat2_index.done"

# use config for tool paths, with defaults
TOOL_CFG = config.get("tools", {}) or {}
STAR_DIR = os.path.expanduser(TOOL_CFG.get("star_dir", "~/pl0296-02/project_data/STAR-2.7.11b/source"))
HISAT2_DIR = os.path.expanduser(TOOL_CFG.get("hisat2_dir", "~/pl0296-02/project_data/hisat2-2.2.1"))
SAMTOOLS_DIR = os.path.expanduser(TOOL_CFG.get("samtools_dir", "~/pl0296-02/project_data/samtools-1.19.2"))
STRINGTIE_DIR = os.path.expanduser(TOOL_CFG.get("stringtie_dir", "~/pl0296-02/project_data/stringtie"))
PYTHON_BIN = os.path.expanduser(TOOL_CFG.get("python_bin", "python"))
SUPPA_BIN = os.path.expanduser(TOOL_CFG.get("suppa2_bin", TOOL_CFG.get("suppa_bin", "suppa.py")))

# default mapping parameters with defaults that can be overridden by sample config
TUNED_PARAMS = config.get("tuned_parameters_defaults", {}) or {}
DEFAULT_MAPPER_PARAMS = {
    "STAR": TUNED_PARAMS.get("star", {}) or {},
    "HISAT2": TUNED_PARAMS.get("hisat2", {}) or {},
}

# Resource defaults for local scheduling inside a single Slurm allocation.
RESOURCE_CFG = config.get("resources", {}) or {}

def resource_mb(key: str, default: int) -> int:
    return int(RESOURCE_CFG.get(key, default) or default)

STAR_MAP_MEM_MB = resource_mb("star_map_mem_mb", 48000)
HISAT2_MAP_MEM_MB = resource_mb("hisat2_map_mem_mb", 12000)
STAR_INDEX_MEM_MB = resource_mb("star_index_mem_mb", 48000)
HISAT2_INDEX_MEM_MB = resource_mb("hisat2_index_mem_mb", 16000)

# Sample design
DATASETS = config.get("datasets", {}) or {}
SAMPLE_ROWS, SAMPLE_IDS = workflowSampleDesign.expand_samples(DATASETS)
            
# Accessors
def sample_cfg(sample: str) -> dict:
    return SAMPLE_ROWS[sample]

def sample_mapper(sample: str) -> str:
    return sample_cfg(sample).get("mapper", "").upper()

def sample_threads(sample: str) -> int:
    return int(sample_cfg(sample).get("threads", 8) or 8)

def sample_mapper_params(sample: str) -> dict:
    mapper = sample_mapper(sample)
    base = DEFAULT_MAPPER_PARAMS.get(mapper, {}).copy() or {}
    base.update(sample_cfg(sample).get("parameters", {}) or {})
    return base

def sample_star_cfg(sample: str) -> dict:
    return sample_mapper_params(sample) if sample_mapper(sample) == "STAR" else {}

def sample_hisat_cfg(sample: str) -> dict:
    return sample_mapper_params(sample) if sample_mapper(sample) == "HISAT2" else {}

def mapper_param(sample: str, key: str, fallback):
    return sample_mapper_params(sample).get(key, fallback)

def star_mm_arg(sample: str) -> str:
    return str(mapper_param(sample, "--outFilterMultimapNmax", 10))

def star_intron_min_arg(sample: str) -> str:
    return str(mapper_param(sample, "--alignIntronMin", 20))

def hisat2_k_arg(sample: str) -> str:
    return str(mapper_param(sample, "-k", 5))

def hisat2_min_intronlen_arg(sample: str) -> str:
    return str(mapper_param(sample, "--min-intronlen", 20))

def sample_param_group(sample: str) -> str:
    return workflowSampleDesign.sample_param_group(sample_cfg(sample))

def mapper_subdir(sample: str) -> str:
    return f"{sample_mapper(sample).lower()}/{sample_param_group(sample)}"

def mapper_bam_pattern(mapper: str) -> str:
    return f"{OUTDIR}/mapping/{mapper.lower()}/{{run}}/{{sample}}.bam"

def mapper_bai_pattern(mapper: str) -> str:
    return f"{mapper_bam_pattern(mapper)}.bai"

def bam_path(sample: str) -> str:
    return f"{OUTDIR}/mapping/{mapper_subdir(sample)}/{sample}.bam"

def sample_output_dir(sample: str) -> str:
    return os.path.dirname(bam_path(sample))

def bai_path(sample: str) -> str:
    return f"{bam_path(sample)}.bai"

def star_read_command(sample: str) -> str:
    explicit = sample_star_cfg(sample).get("readFilesCommand")
    if explicit:
        return explicit

    read_paths = [
        sample_cfg(sample).get("read1", ""),
        sample_cfg(sample).get("read2", ""),
    ]
    compressed_suffixes = (".gz", ".gzip")
    if any(path.endswith(compressed_suffixes) for path in read_paths if path):
        return "zcat"
    return "cat"

def hisat2_read_format_flag(sample: str) -> str:
    read_paths = [
        sample_cfg(sample).get("read1", ""),
        sample_cfg(sample).get("read2", ""),
    ]
    fasta_suffixes = (
        ".fa", ".fa.gz",
        ".fasta", ".fasta.gz",
        ".fna", ".fna.gz",
        ".ffn", ".ffn.gz",
    )
    if any(path.endswith(fasta_suffixes) for path in read_paths if path):
        return "-f"
    return "-q"
STAR_SAMPLES = [s for s in SAMPLE_IDS if sample_mapper(s) == "STAR"]
HISAT2_SAMPLES = [s for s in SAMPLE_IDS if sample_mapper(s) == "HISAT2"]


# Ground truth helpers
def sample_truth_read1(sample: str) -> str:
    return sample_cfg(sample).get("read1")

def sample_truth_read2(sample: str) -> str:
    return sample_cfg(sample).get("read2")

COORDINATE_TARGETS = [
    f"{OUTDIR}/ground_truth/{mapper_subdir(s)}/{s}.coordinates.tsv"
    for s in SAMPLE_IDS
]

STANDARD_SUMMARY_TARGETS = [
    f"{OUTDIR}/ground_truth/{mapper_subdir(s)}/{s}.standard_summary.tsv"
    for s in SAMPLE_IDS
]

GROUND_TRUTH_SUMMARY_TARGETS = [
    f"{OUTDIR}/ground_truth/{mapper_subdir(s)}/{s}.ground_truth_summary.tsv"
    for s in SAMPLE_IDS
]

STRATIFIED_SUMMARY_TARGETS = [
    f"{OUTDIR}/ground_truth/{mapper_subdir(s)}/{s}.stratified_summary.tsv"
    for s in SAMPLE_IDS
]

#ALL_COORDINATES_TARGET = f"{OUTDIR}/ground_truth/all_coordinates.tsv"
ALL_STANDARD_SUMMARY_TARGET = f"{OUTDIR}/ground_truth/all_standard_summary.tsv"
ALL_GROUND_TRUTH_SUMMARY_TARGET = f"{OUTDIR}/ground_truth/all_ground_truth_summary.tsv"
ALL_STRATIFIED_SUMMARY_TARGET = f"{OUTDIR}/ground_truth/all_stratified_summary.tsv"
GROUND_TRUTH_GTF_TABLE_TARGET = f"{OUTDIR}/ground_truth/gtf_exons.tsv"

ANALYSIS_DONE = f"{OUTDIR}/analysis/.done"

MAPPER_QC_TARGETS = [
    f"{OUTDIR}/mapper_qc/{sample_mapper(s).lower()}/{sample_param_group(s)}/{s}.mapper_qc.tsv"
    for s in SAMPLE_IDS
]
ALL_MAPPER_QC_TARGET = f"{OUTDIR}/mapper_qc/all_mapper_qc.tsv"

# AS-event helpers for SUPPA2 annotation on datasets with two transcripts per gene.
AS_EVENT_CFG = config.get("as_events", {}) or {}
AS_EVENT_DATASETS = AS_EVENT_CFG.get("datasets", ["dataset2", "dataset3"]) or []

def dataset_short_name(dataset: str) -> str:
    return dataset.replace("sim_mouse_", "")

def as_event_fasta(dataset: str) -> str:
    return f"data/input/selected_transcripts/transcripts_{dataset_short_name(dataset).replace('dataset', '')}tx.fa"

def as_event_table_target(dataset: str) -> str:
    return f"{OUTDIR}/as_events/{dataset_short_name(dataset)}/as_event_table.tsv"

def as_event_table_targets():
    return [as_event_table_target(dataset) for dataset in AS_EVENT_DATASETS]

def as_event_per_event_table_target(dataset: str) -> str:
    return f"{OUTDIR}/as_events/{dataset_short_name(dataset)}/as_event_table_per_event.tsv"

def as_event_per_event_table_targets():
    return [as_event_per_event_table_target(dataset) for dataset in AS_EVENT_DATASETS]

# rMATS helpers for differential alternative-splicing analysis on datasets 2/3.
RMATS_CFG = config.get("rmats", {}) or {}
RMATS_ENABLED = str(RMATS_CFG.get("enabled", True)).lower() in {"1", "true", "yes", "y"}
RMATS_DATASETS = set(RMATS_CFG.get("datasets", ["sim_mouse_dataset2", "sim_mouse_dataset3"]) or [])
RMATS_CONDITION_A = RMATS_CFG.get("condition_a", {"name": "Cond1", "sample_ids": [1, 2, 3, 4, 5]}) or {}
RMATS_CONDITION_B = RMATS_CFG.get("condition_b", {"name": "Cond2", "sample_ids": [6, 7, 8, 9, 10]}) or {}
RMATS_TRUTH_TABLE_TARGET = RMATS_CFG.get("truth_table") or f"{OUTDIR}/rmats/simulation_truth.tsv"

def sample_number(sample: str) -> int:
    return int(str(sample).rsplit("_", 1)[-1])

def rmats_cases():
    if not RMATS_ENABLED:
        return []
    cases = set()
    for sample in SAMPLE_IDS:
        row = sample_cfg(sample)
        dataset = row.get("dataset", "")
        if dataset in RMATS_DATASETS:
            cases.add((dataset, sample_mapper(sample).lower(), sample_param_group(sample)))
    return sorted(cases)

def rmats_bams(dataset: str, mapper: str, param_group: str, condition: dict) -> list:
    wanted_ids = {int(x) for x in condition.get("sample_ids", [])}
    samples = [
        sample for sample in SAMPLE_IDS
        if sample_cfg(sample).get("dataset", "") == dataset
        and sample_mapper(sample).lower() == mapper.lower()
        and sample_param_group(sample) == param_group
        and sample_number(sample) in wanted_ids
    ]
    return [bam_path(sample) for sample in sorted(samples, key=sample_number)]

def rmats_summary_targets():
    return [
        f"{OUTDIR}/rmats/{dataset}/{mapper}/{param_group}/summary.tsv"
        for dataset, mapper, param_group in rmats_cases()
    ]

def rmats_all_summary_target():
    return f"{OUTDIR}/rmats/all_summary.tsv"

def rmats_truth_design_files():
    dataset_names = sorted(RMATS_DATASETS)
    return [
        f"data/input/sim/polyester_design/{dataset_short_name(dataset)}/transcript_design.tsv"
        for dataset in dataset_names
    ]

def rmats_targets():
    targets = rmats_summary_targets()
    if targets:
        targets.append(rmats_all_summary_target())
    return targets


def star_bam_targets():
    return [bam_path(s) for s in STAR_SAMPLES]

def star_bai_targets():
    return [bai_path(s) for s in STAR_SAMPLES]

def hisat2_bam_targets():
    return [bam_path(s) for s in HISAT2_SAMPLES]

def hisat2_bai_targets():
    return [bai_path(s) for s in HISAT2_SAMPLES]

def sample_design_target():
    return f"{OUTDIR}/sample_design.tsv"

include: "rules/refs.smk"
include: "rules/mapping.smk"
include: "rules/sample_design.smk"
include: "rules/as_events.smk"
include: "rules/ground_truth.smk"
include: "rules/rmats.smk"
include: "rules/mapper_qc.smk"
include: "rules/analysis.smk"

rule all:
    input:
        STAR_INDEX_DONE,
        HISAT2_INDEX_DONE,
        star_bam_targets(),
        star_bai_targets(),
        hisat2_bam_targets(),
        hisat2_bai_targets(),
        sample_design_target(),
        as_event_table_targets(),
        as_event_per_event_table_targets(),
        GROUND_TRUTH_GTF_TABLE_TARGET,
        COORDINATE_TARGETS,
        STANDARD_SUMMARY_TARGETS,
        GROUND_TRUTH_SUMMARY_TARGETS,
        STRATIFIED_SUMMARY_TARGETS,
        #ALL_COORDINATES_TARGET,
        ALL_STANDARD_SUMMARY_TARGET,
        ALL_GROUND_TRUTH_SUMMARY_TARGET,
        ALL_STRATIFIED_SUMMARY_TARGET,
        rmats_targets(),
        MAPPER_QC_TARGETS,
        ALL_MAPPER_QC_TARGET,
        ANALYSIS_DONE,
