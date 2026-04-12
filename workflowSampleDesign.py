import os
import re

SAMPLE_DESIGN_HEADERS = [
    "sample",
    "dataset",
    "param_group",
    "mapper",
    "read1",
    "read2",
    "threads",
    "strandedness",
]

def sample_param_group(row: dict) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(row.get("param_group", "")).strip())
    safe = re.sub(r"_+", "_", safe).strip("_")
    return safe or "run"

def expand_samples(datasets: dict) -> tuple[dict, list]:
    sample_rows = {}
    sample_ids = []

    for dataset, ds_cfg in (datasets or {}).items():
        read1_tmpl = ds_cfg.get("read1", "") or ""
        read2_tmpl = ds_cfg.get("read2", "") or ""
        dataset_sample_ids = ds_cfg.get("sample_ids", []) or []
        param_groups = ds_cfg.get("param_groups", {}) or {}
        dataset_threads = ds_cfg.get("threads", 8) or 8

        for group_name, group_cfg in param_groups.items():
            mapper = str(group_cfg.get("mapper", "") or "").upper()
            params = group_cfg.get("parameters", {}) or {}
            threads = group_cfg.get("threads", dataset_threads)
            strandedness = group_cfg.get("strandedness", ds_cfg.get("strandedness", "") or "")

            for sample_id in dataset_sample_ids:
                sample = f"{group_name}_{sample_id}"
                read1 = read1_tmpl.format(sample_id=sample_id)
                read2 = read2_tmpl.format(sample_id=sample_id) if read2_tmpl else ""
                sample_rows[sample] = {
                    "sample": sample,
                    "dataset": dataset,
                    "param_group": group_name,
                    "mapper": mapper,
                    "read1": read1,
                    "read2": read2,
                    "threads": threads,
                    "strandedness": strandedness,
                    "parameters": params,
                }
                sample_ids.append(sample)

    return sample_rows, sample_ids


def write_sample_design(path: str, sample_rows: dict, sample_ids: list) -> None:
    outdir = os.path.dirname(path)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    with open(path, "w") as handle:
        handle.write("\t".join(SAMPLE_DESIGN_HEADERS) + "\n")
        for sample in sorted(sample_ids):
            row = sample_rows[sample]
            values = [
                str(row.get("sample", "")),
                str(row.get("dataset", "")),
                str(row.get("param_group", "")),
                str(row.get("mapper", "")),
                str(row.get("read1", "")),
                str(row.get("read2", "")),
                str(row.get("threads", "")),
                str(row.get("strandedness", "")),
            ]
            handle.write("\t".join(values) + "\n")
