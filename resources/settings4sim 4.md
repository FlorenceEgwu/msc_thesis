# Settings for datasets generation

## 1 Syntax of `simulate_experiment`

We are using `simulate_experiment` function with syntax:

```r
library(polyester)
simulate_experiment(fasta_path, reads_per_transcript = readspertx,
    num_reps = c(5, 5), fold_changes = fold_changes, outdir = "choose_your_place")
```

where

- `fasta_path` is a path to `cdna.fa` file in which are stored sequences of chosen transcripts,
- `readspertx` has the following form:

```r
readspertx = round(cov_lvl * width(fasta)/read_length)
```

We set `read_length = 100`, where `fasta` is an object loaded from `fasta_path` with *Biostrings* package and `cov_lvl` is a vector of values `5, 10, 20` that will be described below in sec. Section 2.2 and Section 3.2

- `fold_changes` we also describe below in sec. Section 2.2 and Section 3.2

## 2 **Dataset 1**

### 2.1 Choosing reference transcripts

We have *900* non-overlapping genes with *1 transcript* with *300 tr with =< 2 exons* and *600 tr with > 2 exons.*

The transcripts are ordered in the form:

1. Tr with =< 2 exons (easy tr)…
  300. Tr with =< 2 exons (easy tr)
  301. Tr with > 2 exons (complex tr)…
    900. Tr with > 2 exons (complex tr)

The structure is `easy_tr x 300 + complex_tr x 600`, so in total we have 900 tr in order that will help during simulation.

Sequences of chosen transcripts are stored in object `fasta_file`of class `XStringSet`. This object can be obtained from fasta file of *cdna* downloaded from ensembl.

```r
library(Biostrings)
fasta_file = "my_sequnces.cdna.fa"
fasta = readDNAStringSet(fasta_file)
```

If you have a list of names of th transcripts you want to use, you can load full `.cda.fa` file and filter it by names of transcripts. Access to names of `XStringSet` object you get with `names(fasta)`.

### 2.2 Simulations settings

For *Dataset 1* ordered in this way we still need to specify values of `reads_per_transcript` and `fold_changes` parameters. Since chosen the dataset is divided in 300 easy and 600 compex transcripts the number of repetitions is changed to reflect that that:

| type of transcript | `fold_changes` | `cov_lvl` | number of repetitions |
| --- | --- | --- | --- |
| easy | 1:1 | 5 | x40 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 1:2 | 5 | x15 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 2:1 | 5 | x15 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 1:4 | 5 | x15 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 4:1 | 5 | x15 |
|  |  | 10 |  |
|  |  | 20 |  |
| complex | 1:1 | 5 | x80 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 1:2 | 5 | x30 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 2:1 | 5 | x30 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 1:4 | 5 | x30 |
|  |  | 10 |  |
|  |  | 20 |  |
|  | 4:1 | 5 | x30 |
|  |  | 10 |  |
|  |  | 20 |  |

To check if the table fir the 900 transcripts in Dataset1 we calculate number of cases. For first 300 transcripts with easy structure we have:

\(3\cdot 40 \;(\rm{tr\;with\;1:1\;fc}) + 3 \cdot 15\cdot 2 \;(\rm{tr\;with\;1:2\;or\;2:1\;fc}) + 3\cdot 15\cdot 2 \;(\rm{tr\;with\;1:4\;or\;4:1\;fc}) =300\;tr\)

For remaining 600 transcripts with complex structure we have:

\(3\cdot 80 \;(\rm{tr\;with\;1:1\;fc}) + 3\cdot 30\cdot 2 \;(\rm{tr\;with\;1:2\;or\;2:1\;fc}) + 3\cdot 30\cdot 2 \;(\rm{tr\;with\;1:4\;or\;4:1\;fc}) =600\;tr\)

## 3 **Dataset 2**

### 3.1 Choosing reference transcripts

We have *900* non-overlapping genes with *2 transcripts* (*1800* total tr) with *300 genes with =< 2 exons* and *600 genes with > 2 exons.*

From the photo you send me all the first tr in genes are supposed to have fold change 1-1 and for the second transcripts you have exactly the same structure in simulation like for Dataset 1. That is why it is usuful to order then as follows:

1. First tr in gene 1 with =< 2 exons (easy tr)
2. First tr in gene 2 with =< 2 exons (easy tr)…
  300. First tr in gene 300 with =< 2 exons (easy tr)
  301. First tr in gene 301 with > 2 exons (complex tr)…
    900. First tr in gene 900 with > 2 exons (complex tr)
3. Second tr in gene 1 with =< 2 exons (easy tr)
4. Second tr in gene 2 with =< 2 exons (easy tr)

```
…

1200. Second tr in gene 300 with =\< 2 exons (easy tr)

1201. Second tr in gene 301 with \> 2 exons (complex tr)

     …

     1800. Second tr in gene 900 with \> 2 exons (complex tr)
```

### 3.2 Simulations settings

For *Dataset 2* first 900 tr will have foldchange 1-1 so we only need to include levels of coverage and division on easy/complex genes. The second part will remain exactly the same like for the previous case:

| which tr in pair | type of transcript | `fold_changes` | `cov_lvl` | number of repetitions |
| --- | --- | --- | --- | --- |
| 1 | easy | 1:1 | 5 | x100 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  | complex | 1:1 | 5 | x200 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
| 2 | easy | 1:1 | 5 | x40 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:2 | 5 | x15 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:1 | 5 | x15 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:4 | 5 | x15 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:1 | 5 | x15 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  | complex | 1:1 | 5 | x80 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:2 | 5 | x30 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:1 | 5 | x30 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:4 | 5 | x30 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:1 | 5 | x30 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |

To check if the table for the second in pair contain 900 transcripts in Dataset2 we calculate number of cases. For first 300 transcripts with easy structure we have:

\(3\cdot 40 \;(\rm{tr\;with\;1:1\;fc}) + 3\cdot 15 \cdot 2 \;(\rm{tr\;with\;1:2\;or\;2:1\;fc}) + 3 \cdot 15 \cdot 2 \;(\rm{tr\;with\;1:4\;or\;4:1\;fc}) =300\;tr\)

For remaining 600 transcripts with complex structure we have:

\(3\cdot 80 \;(\rm{tr\;with\;1:1\;fc}) + 3 \cdot 30 \cdot 2 \;(\rm{tr\;with\;1:2\;or\;2:1\;fc}) + 3\cdot 30 \cdot 2 \;(\rm{tr\;with\;1:4\;or\;4:1\;fc}) =600\;tr\)

## 4 **Dataset 3**

### 4.1 Choosing reference transcripts

We again have *900* non-overlapping genes with *2 transcripts* (*1800* total tr) with *300 genes with =< 2 exons* and *600 genes with > 2 exons.*

Lets organize transcripts in the same order as for dataset 2:

1. First tr in gene 1 with =< 2 exons (easy tr)
2. First tr in gene 2 with =< 2 exons (easy tr)…
  300. First tr in gene 300 with =< 2 exons (easy tr)
  301. First tr in gene 301 with > 2 exons (complex tr)…
    900. First tr in gene 900 with > 2 exons (complex tr)
3. Second tr in gene 1 with =< 2 exons (easy tr)
4. Second tr in gene 2 with =< 2 exons (easy tr)

```
…

1200. Second tr in gene 300 with =\< 2 exons (easy tr)

1201. Second tr in gene 301 with \> 2 exons (complex tr)

     …

     1800. Second tr in gene 900 with \> 2 exons (complex tr)
```

### 4.2 Simulations settings

For *Dataset 3* first 900 tr (first tr in pair) will have foldchanges with increasing expression in second group (ie 1-1, 1-2, 1-4, 2-4) . The second part (second tr in pair) will have foldchanges with increasing expression in first group (ie 1-1, 2-1, 4-1, 4-2):

| which tr in pair | type of transcript | `fold_changes` | `cov_lvl` | number of repetitions |
| --- | --- | --- | --- | --- |
| 1 | easy | 1:1 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:2 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:4 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:4 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  | complex | 1:1 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:2 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 1:4 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:4 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
| 2 | easy | 1:1 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:1 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:1 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:2 | 5 | x25 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  | complex | 1:1 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 2:1 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:1 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |
|  |  | 4:2 | 5 | x50 |
|  |  |  | 10 |  |
|  |  |  | 20 |  |

To check if the table contain proper number of transcripts (900 for first in pair and 900 for second in pair) we calculate number of cases for the first in pair transcripts. Second in pair have exactly the same structure. For first 300 transcripts with easy structure we have:

\(3 (\#\;of \;cases \;for \;cov)\cdot 25 (\#\;of \;repetitions)\cdot 4 (\#\;of \;cases \; for \;fc) =300\;tr\)

For remaining 600 transcripts with complex structure we have:

\(3 (\#\;of \;cases \;for \;cov)\cdot 50 (\#\;of \;repetitions)\cdot 4 (\#\;of \;cases \; for \;fc) =600\;tr\)

## 5 Mapping settings

We are focus on 2 parameters with 3 possible options.

In STAR they are:

- `outFilterMultimapNmax` with values 5, 10, 20
- `alignIntronMin` with values 20 50 100

In Hisat2 they are:

- `-k`
- `--min-intronlen`

with the same values as for STAR.

## 6 Standard summary of mapping

Standard summuaries of mapping (eg number of reads in analysis, number of mapped reads, % uniquely mapped reads, % reads mapped to multiple loci, % unmapped reads etc) can be obtained from files generated by mappers. For STAR they are stored in files `Log.final.out`. For hisat2 it is stored in log files.

## 7 Finding “ground truth” mapping statistics

We will use several tables to get this statistics.

1. Firstly, you need Table GTF of all trancripts that were used in the analysis (from gtf) that will include start, end, strand of each exon in transcripts.
2. Secondly, you need Table COORDINATES based on bam file with coordinates for each read, such as chromosome, start, end, strand, cigar (they are stored columns called RNAME, POS, CIGAR) and the real coordinates of this read (stored in QNAME).

- For each bam file you will have separate table.
- The information in QNAME includes real positions of read based on transcriptome (not genome) so it needs to be conversed to genome position (using gtf), so mapped and real locations could be compared.

3. Next, based on Table COORDINATES you should check which reads were alinged in correct transcripts (using for example GRange or GenomicRanges with functions such as findOverlaps). You can include this result as separate column in the Table COORDINATES as well as additional column with the information based on CIGAR if the read was spliced (includes juction) and the information from polyester simulation if this transcript was easy or complex (number of exon is <= 3).
4. Finally, you calculate the statistics in Table MAPPING GROUND TRUTH for each sample in each considered case 540 rows (3 datasets x 10 samples in datasets x 2 mappers x 3 levels of mutlimap parameter x 3 levels of intron size parameter). In this table you can have columns describing the accuracy of mapping. In addition you can also ganerate additional tables with the same summaries but with the division based on:

- type of read (spliced/simple);
- type of mapping error (eg incorrect_position - read in completely wrong region vs incorrect_junction - read in transcript location, but wrong junctions or located outside exons);
- type of transcript structue (easy/complex transcipt);

The results of Tables MAPPING GROUND TRUTH will be starting points for visualizations and tables that should be included in the Results section of the thesis.

### 7.1 Division of summarizations

We would like to see how the mappers with different parameters behave in different cases. We can include here several cases:

- type of transcript structue: Do we have higher precision for “easy” genes compared to “complex” ones?
- type of read: Do we have higher precision for “simple” reads compared to “spliced” ones?
- type of mapping error: How many mismapped reads (percentages) are in completely wrong region and how many of them are close proper loaction. I would generate here the results as a table in which you would present % of mismapped read for mapping error
- type of alternative splicing event. Summarisations here can be done for dataset 2 and 3 only. For this division you need to first apply SUPPA2 or IsoformSwitchAnalyzeR package to get the information for each gene (with 2 transcripts) chosen for simulation what alternative splicing event it represent. This information should be included in the table you ganerate for summarizations. For this part we can draw plot with division on AS type to see if any mapper or any parameter’s settings result in lower precision for some particular AS events.

## 8 rMATS analysis

The last type of the analysis would be checking changes in differential alvernative splicing for different mappers and its parameter’s settings for datasets 2 and 3.

For each case that we have (2 dataset x 2 mappers x 3 levels of mutlimap parameter x 3 levels of intron size parameter) 36 cases rMATS needs to be applied. For each case list of isoforms with significant changes would be generated. We would need to check if the found isoforms agrees with setting of fc in simulations.
