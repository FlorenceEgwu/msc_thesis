# ---------------------------
# Installation of dependencies
# ---------------------------
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

# Ensure BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
install_if_missing("devtools")
if (!requireNamespace("polyester", quietly = TRUE)) {
  devtools::install_github("alyssafrazee/polyester", dependencies = TRUE)
}
install_if_missing("Biostrings", bioc = TRUE)
install_if_missing("R.utils")

library(polyester)
library(Biostrings)
library(R.utils)


# ---------------------------
# User-adjustable parameters
# ---------------------------
fold_changes = matrix(c(4,4,rep(1,18),1,1,4,4,rep(1,16)), nrow=20)

# FASTA annotation
fasta_file = system.file('extdata', 'chr22.fa', package='polyester')
fasta = readDNAStringSet(fasta_file)

# subset the FASTA file to first 20 transcripts
small_fasta = fasta[1:20]
writeXStringSet(small_fasta, 'chr22_small.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM
readspertx = round(20 * width(small_fasta) / 100)

# simulation call: 
simulate_experiment('chr22_small.fa', reads_per_transcript=readspertx, 
                    num_reps=c(10,10), fold_changes=fold_changes, outdir='simulated_reads')
