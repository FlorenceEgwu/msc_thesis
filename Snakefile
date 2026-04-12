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

rule all:
    input:
        STAR_INDEX_DONE,
        HISAT2_INDEX_DONE,
        star_bam_targets(),
        star_bai_targets(),
        hisat2_bam_targets(),
        hisat2_bai_targets(),
        sample_design_target(),
