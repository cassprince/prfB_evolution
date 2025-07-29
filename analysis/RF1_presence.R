library(jsonlite)
library(tidyverse)
library(readr)

setwd("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//prfB//Data//Bioinformatics")

df_GCF_nc = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//GCF_nuccore_reps_clean.csv")

lines = readLines("assembly_data_report.jsonl")
lines = lapply(lines, fromJSON)
lines = lapply(lines, unlist)
df_GCF_tax = bind_rows(lines)

df_GCF_tax_select = df_GCF_tax %>%
  select(accession, assemblyInfo.assemblyLevel, organism.organismName, assemblyInfo.bioprojectLineage.bioprojects.title, organism.taxId, checkmInfo.checkmMarkerSet, checkmInfo.completeness, checkmInfo.contamination, assemblyStats.totalSequenceLength, assemblyStats.gcPercent, assemblyStats.numberOfComponentSequences)

#############
# Join the NCBI accessions with their metadata.
df_GCF_nc_tax = left_join(df_GCF_nc, df_GCF_tax_select, by = join_by(assembly == accession)) %>%
  mutate(checkmInfo.completeness = as.numeric(checkmInfo.completeness)) %>%
  mutate(checkmInfo.contamination = as.numeric(checkmInfo.contamination)) %>%
  mutate(assemblyStats.totalSequenceLength = as.numeric(assemblyStats.totalSequenceLength)) %>%
  mutate(assemblyStats.gcPercent = as.numeric(assemblyStats.gcPercent)) %>%
  mutate(assemblyStats.numberOfComponentSequences = as.numeric(assemblyStats.numberOfComponentSequences))
##############

lineages = read.csv("C://Users//cassp//Box Sync//Feaga Lab//Cassidy Prince//Katrina//ref_lineage.txt", col.names = "taxID", header = FALSE) %>% 
  separate(taxID, into = c("organism", "domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";", extra = "merge") %>% 
  separate(organism, into = c("taxID", "organism"), sep = "\\s", extra = "merge") %>%
  drop_na()


df = read.csv("dataset.txt", header = FALSE, col.names = "ID") 
df_tax = inner_join(df_GCF_tax_select, lineages, by = join_by("organism.taxId" == "taxID")) %>%
  filter(accession %in% df$ID) %>%
  distinct()


## Clean prfA 

df_prfA = read.csv("prfA_desc.txt", header = FALSE, col.names = "description") 

df_prfA = df_prfA %>%
  mutate(ID = substr(description, 6, nchar(description))) %>%
  mutate(prfA_pseudo = grepl("pseudo", df_prfA$description))

df_prfA = df_prfA %>%
  mutate(ID, ID = sub('\\_cds.*', '', df_prfA$ID))

df_prfA_GCF = inner_join(df_prfA, df_GCF_nc, by = join_by("ID" == "nuccore")) 

## Clean RF1

df_RF1 = read.csv("RF1_desc.txt", header = FALSE, col.names = "description") 

df_RF1 = df_RF1 %>%
  mutate(ID = substr(description, 6, nchar(description))) %>%
  mutate(RF1_pseudo = grepl("pseudo", df_RF1$description))

df_RF1 = df_RF1 %>%
  mutate(ID, ID = sub('\\_cds.*', '', df_RF1$ID))

df_RF1_GCF = inner_join(df_RF1, df_GCF_nc, by = join_by("ID" == "nuccore")) 

## Clean RF_like

df_RF_like = read.csv("RF_like_desc.txt", header = FALSE, col.names = "description") 

df_RF_like = df_RF_like %>%
  mutate(ID = substr(description, 6, nchar(description))) %>%
  mutate(RF_like_pseudo = grepl("pseudo", df_RF_like$description))

df_RF_like = df_RF_like %>%
  mutate(ID, ID = sub('\\_cds.*', '', df_RF_like$ID))

df_RF_like_GCF = inner_join(df_RF_like, df_GCF_nc, by = join_by("ID" == "nuccore")) 


### Join dataset with RF data to make final dataset

df_final = df_tax %>%
  mutate(prfA = df_tax$accession %in% df_prfA_GCF$assembly) %>% 
  mutate(RF1 = df_tax$accession %in% df_RF1_GCF$assembly) %>% 
  mutate(RF_like = df_tax$accession %in% df_RF_like_GCF$assembly)

df_final_pseudos = df_final %>% 
  left_join(select(df_prfA_GCF, assembly, prfA_pseudo), by = join_by("accession" == "assembly")) %>% 
  left_join(select(df_RF1_GCF, assembly, RF1_pseudo), by = join_by("accession" == "assembly")) %>% 
  left_join(select(df_RF_like_GCF, assembly, RF_like_pseudo), by = join_by("accession" == "assembly")) 

### Summarize

yes_RF1 = df_final %>% filter(prfA == TRUE | RF1 == TRUE) %>% distinct(accession, .keep_all = TRUE)
RF1_like = df_final %>% filter(RF_like == TRUE) %>% distinct(accession, .keep_all = TRUE)
RF1_like_only = df_final %>% filter(prfA == FALSE & RF1 == FALSE & RF_like == TRUE) %>% distinct(accession, .keep_all = TRUE)
no_RF1 = df_final %>% filter(prfA == FALSE & RF1 == FALSE & RF_like == FALSE) %>% distinct(accession, .keep_all = TRUE)

RF1_like_only = df_final %>% filter(prfA == FALSE & RF1 == FALSE & RF_like == TRUE)

# Are there any versions of RF1 that are pseudogenized
pseudos = df_final_pseudos %>% filter(prfA_pseudo == TRUE | RF1_pseudo == TRUE | RF_like_pseudo == TRUE) 


pseudos_one = pseudos %>% group_by(accession) %>% filter(n() == 1)

RF1_pseudos = pseudos_one %>% 
  filter(prfA == TRUE & RF1 == TRUE) %>% 
  filter(prfA_pseudo == TRUE | RF1_pseudo == TRUE) 


prfA_pseudos = pseudos %>% filter(RF1_pseudo == FALSE)

