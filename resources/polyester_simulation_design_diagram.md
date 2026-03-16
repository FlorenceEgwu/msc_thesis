# Polyester Simulation Design Diagram

Run output base:
`results/polyester_design_run_2026-02-27_final`

```mermaid
flowchart TD
  A[Input FASTAs] --> B1[dataset_1: 900 genes, 900 transcripts]
  A --> B2[dataset_2: 900 genes, 1800 transcripts]

  B1 --> C1[rpt split: 5/10/20 = 300/300/300 genes]
  B1 --> C2[fc split: 1/1, 1/2, 1/4 = 300/300/300 genes]
  B1 --> C3[direction: equal/A_gt_B/B_gt_A = 300/300/300]
  B1 --> O1[Outputs: 4 paired FASTQ + gene_design.tsv + transcript_design.tsv + sample_design.tsv + count_matrix.tsv]

  B2 --> D1[rpt split: 5/10/20 = 300/300/300 genes]
  B2 --> D2[fc split: 1/1, 1/2, 1/4 = 300/300/300 genes]
  B2 --> D3[direction: equal/A_gt_B/B_gt_A = 300/300/300]
  B2 --> O2[Outputs: 4 paired FASTQ + gene_design.tsv + transcript_design.tsv + sample_design.tsv + count_matrix.tsv]
```

## Plain Text Version

- dataset_1
  - genes: 900
  - transcripts: 900
  - rpt levels: 5, 10, 20 (300 genes each)
  - fc ratios: 1:1, 1:2, 1:4 (300 genes each)
  - directions: equal, A_gt_B, B_gt_A (300 genes each)
  - outputs: sample_01_1/2.fasta.gz, sample_02_1/2.fasta.gz, gene_design.tsv, transcript_design.tsv, sample_design.tsv, count_matrix.tsv

- dataset_2
  - genes: 900
  - transcripts: 1800
  - rpt levels: 5, 10, 20 (300 genes each)
  - fc ratios: 1:1, 1:2, 1:4 (300 genes each)
  - directions: equal, A_gt_B, B_gt_A (300 genes each)
  - outputs: sample_01_1/2.fasta.gz, sample_02_1/2.fasta.gz, gene_design.tsv, transcript_design.tsv, sample_design.tsv, count_matrix.tsv
