# Ground-Truth Mapping Notes

This module is designed to work without changing any existing workflow files.

It writes new outputs under `data/results/ground_truth/`:

- `gtf_transcript_exons.tsv.gz`
- `read_truth/dataset*_sample_*.tsv.gz`
- `coordinates/<sample>.tsv.gz`
- `mapping_ground_truth.tsv.gz`
- grouped summary tables

## Important assumption

The current BAM `QNAME` keeps only the `readNNN` prefix, while the original Polyester FASTA headers contain the transcript and transcript-space coordinates.

Example:

- FASTA header: `read297420/ENSMUST...;mate1:1-100;mate2:157-255`
- BAM `QNAME`: `read297420`

Because the mate-level coordinate suffix was dropped before alignment, the script:

1. joins BAM records back to the original Polyester header by `read_id`
2. lifts both header coordinate candidates (`mate1`, `mate2`) from transcript space to genome space
3. assigns the BAM record to the closer of the two lifted candidates

This is exact for most correctly aligned reads and is the safest approximation available from the current files.

## Transcript structure columns

The coordinates and summary tables include both:

- `exon_structure_type`: derived from the GTF using `exon_count <= exon_threshold` as `easy`
- `sim_transcript_type`: the label recorded in the Polyester transcript design table

This keeps the output usable even if the simulation label and the GTF-based threshold are not identical.
