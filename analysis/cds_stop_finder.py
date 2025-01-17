# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:32:31 2023

@author: cassp
"""

from Bio import SeqIO 
import pandas as pd
import re

file = SeqIO.parse(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\cds_myco.fna", "fasta")

#Determine the number of CDS in the multifasta file.
num = len([1 for line in open(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\cds_myco.fna") if line.startswith(">")])

df = pd.DataFrame(None, index = range(num), columns=["genome_ID", "gene", "protein", "locus_tag", "terminal_stop"])

count = 0
for seq_record in file:
        seq = seq_record.seq
        ID_result = re.search(r"lcl\|(.*?)\_cds", seq_record.id)
        print(ID_result.group(1))
        
        ID = ID_result.group(1)
        sub_ID = ID
        locus_result = re.search(r"locus_tag\=(.*?)\]", seq_record.description)

        df.loc[count, "genome_ID"] = sub_ID
        df.loc[count, "locus_tag"] = locus_result.group(1)

        if "{gene=" in seq_record.description:
            gene_result = re.search(r"gene\=(.*?)\]", seq_record.description)
            df.loc[count, "gene"] = gene_result.group(1)

        if "protein" in seq_record.description:
            protein_result = re.search(r"protein\=(.*?)\]", seq_record.description)
            df.loc[count, "protein"] = protein_result.group(1)

        df.loc[count, "terminal_stop"] = seq[len(seq)-3:]

        count += 1

        
df.to_csv(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\cds_myco_stops_2.csv")






