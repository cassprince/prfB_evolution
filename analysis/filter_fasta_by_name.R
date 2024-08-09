# Load in packages, fasta file with prfB sequences, and dataframe containing only genome accessions with <10% CheckM contamination.
library(tidyverse)
library(seqinr)
library(readxl)

file = "C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/prfB/Data/Bioinformatics/prfB_NCBI_no_dups_RC.fasta"
df = data.frame(read_csv("C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/prfB/Data/Bioinformatics/FS_data_clean_8_1_24.csv"))
accs = df[,"nuccore"]

fasta = read.fasta(file, as.string = TRUE)
annot = getAnnot(fasta)
output = list()

# Filter prfB sequences to make sure only genomes with <10% CheckM contamination are included.
for (n in accs){
  output = append(output, fasta[grep(n, annot)])
}

write.fasta(output, names = names(output), file.out = "C:/Users/cassp/Box Sync/Feaga Lab/Cassidy Prince/prfB/Data/Bioinformatics/prfB_filtered.fasta")
