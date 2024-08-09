# Load packages and set working directory.

library(jsonlite)
library(tidyverse)
library(readr)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//prfB//Data//Bioinformatics")

# Upload prfB frameshift data created from Python script. Remove location of gene hit from ID column to get the NCBI Nucleotide accession number.
df_FS = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//prfB//Data//Bioinformatics//FS_data_8_5_24.csv") 

df_FS = df_FS %>%
  mutate(ID, ID = gsub("[^:]+$", "", df_FS$ID)) 
df_FS = df_FS %>%
  mutate(ID, ID = gsub(":", "", df_FS$ID)) %>%
  select(-X, -prfB_DNA_seq, -prfB_DNA_seq_alignment, -prfB_AA_seq, -seq_description, -notes)

# Upload lineages acquired from NCBI taxdump and taxonkit based on each assembly's taxid.
lineages = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ref_lineage.txt", col.names = "taxID", header = FALSE) %>% 
  separate(taxID, into = c("organism", "domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>% 
  separate(organism, into = c("taxID", "organism"), sep = "\\s", extra = "merge") %>%
  drop_na()

# Upload NCBI RefSeq Assembly database accessions (beginning with "GCF_") and corresponding NCBI Nucleotide (nuccore) database accessions (beginning with "NC_" or "NZ_".
df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")

# Upload assembly metadata for NCBI assemblies and select relevant data columns.
lines = readLines("assembly_data_report.jsonl")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_GCF_tax = bind_rows(lines)

df_GCF_tax_select = df_GCF_tax %>%
  select(accession, assemblyInfo.assemblyLevel, organism.organismName, assemblyInfo.bioprojectLineage.bioprojects.title, organism.taxId, checkmInfo.checkmMarkerSet, checkmInfo.completeness, checkmInfo.contamination, assemblyStats.totalSequenceLength, assemblyStats.gcPercent, assemblyStats.numberOfComponentSequences)

# Join the NCBI accessions with their metadata.
df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == accession)) %>%
  mutate(checkmInfo.completeness = as.numeric(checkmInfo.completeness)) %>%
  mutate(checkmInfo.contamination = as.numeric(checkmInfo.contamination)) %>%
  mutate(assemblyStats.totalSequenceLength = as.numeric(assemblyStats.totalSequenceLength)) %>%
  mutate(assemblyStats.gcPercent = as.numeric(assemblyStats.gcPercent)) %>%
  mutate(assemblyStats.numberOfComponentSequences = as.numeric(assemblyStats.numberOfComponentSequences))

# Prepare the final dataframe with Assembly and Nucleotide accession numbers, assembly statistics, taxonomy, and gene presence.
# Remove archaeal genomes. Filter based on CheckM contamination. Remove any duplicate genomes.
df_distinct = inner_join(df_FS, df_GCF_nc_tax, by = join_by("ID" == "nuccore")) %>%
  select(-X) %>%
  rename(nuccore = ID) %>% 
  inner_join(lineages, by = join_by("organism.taxId" == "taxID"), multiple = "all") %>%
  filter(!grepl("Archaea", domain)) %>%
  filter(checkmInfo.contamination < 10) %>%
  distinct(assembly, .keep_all = TRUE) %>%  
  relocate(assembly, nuccore)

# Rename taxonomy to be consistent with Coleman et al 2021.
df_distinct$phylum[df_distinct$phylum == "delta/epsilon subdivisions"] = "Pseudomonadota"
df_distinct$phylum[df_distinct$phylum == "Pseudomonadota"] = "Proteobacteria"
df_distinct$phylum[df_distinct$phylum == "Terrabacteria group"] = df_distinct$class[df_distinct$phylum == "Terrabacteria group"]
df_distinct$phylum[df_distinct$phylum == "Bacillota"] = "Firmicutes"
df_distinct$phylum[df_distinct$phylum == "Actinomycetota"] = "Actinobacteriota"
df_distinct$phylum[df_distinct$phylum == "Abditibacteriota"] = "Armatimonadota"
df_distinct$phylum[df_distinct$phylum == "Aquificota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Campylobacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Deferribacterota"] = "Aquificota + Campylobacterota + Deferribacterota"
df_distinct$phylum[df_distinct$phylum == "Thermodesulfobacteriota"] = "Desulfuromonadota + Desulfobacterota"
df_distinct$phylum[df_distinct$phylum == "Calditrichota"] = "FCB group"
df_distinct$phylum[df_distinct$phylum == "Thermomicrobiota"] = "Chloroflexota"
df_distinct$phylum[df_distinct$phylum == "Proteobacteria"] = df_distinct$class[df_distinct$phylum == "Proteobacteria"]

# Write dataframe for figures.
write.csv(df_distinct, "FS_data_clean_8_5_24.csv")
