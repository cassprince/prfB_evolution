# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 15:03:00 2023

@author: cassp
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



file = SeqIO.parse(r"C:\Users\cassp\OneDrive\Documents\Feaga Lab\prfB\prfB_aab.fasta", "fasta")
seq_list = []
seq_des_list = []

for seq_record in file:
    seq = seq_record.seq
    print(seq_record)
    #print(seq_record.description)
    if seq[1:3] != "TG":
        rec = SeqRecord(seq.reverse_complement(), id = seq_record.id, description= seq_record.description)
        seq_list.append(rec)
    else:
        rec = SeqRecord(seq, id = seq_record.id, description= seq_record.description)
        seq_list.append(rec)

SeqIO.write(seq_list, r"C:\Users\cassp\OneDrive\Documents\Feaga Lab\prfB\prfB_aab_RC.fasta", "fasta")