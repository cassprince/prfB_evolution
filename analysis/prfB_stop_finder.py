# -*- coding: utf-8 -*-
"""
Created on Sat Mar 11 14:42:29 2023

@author: cassp
"""

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import re

alignment = AlignIO.read(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\prfB_NCBI_filt_align.fasta", "fasta")

df = pd.DataFrame(0, index = range(len(alignment)), columns = ["ID", "start", "in_frame_stop?", "internal_stop_index", "SD_seq", "FS_seq", "internal_stop", "nuc_after_stop", "terminal_stop", "FS_context_seq", "prfB_DNA_seq", "prfB_DNA_seq_alignment", "prfB_AA_seq", "seq_description", "notes"])
df_nd = pd.DataFrame(0, index = range(len(alignment)), columns = ["ID"])

expectedStop = 527
table = 11
stops = ["TAA", "TGA", "TAG"]

count = 0
for seq in alignment:
    print(seq.id)
    #print(seq.description)
    
    seqString = str(seq.seq)
    seqString = seqString.upper()
    seqStringNH = seqString.replace("-", "")
    seqStringToTranslate = Seq(seqStringNH)
    ID = seq.id
    shortID = ID.split(':', 1)[0]
    
#Identify sequences that lack the expected region for internal stop codons.
    if all(ch == "-" for ch in seq[:expectedStop]) or seqStringNH[:3] == "---":         ##If everything before the expected location of the TGA based on E. coli prfB is hyphens, record that the sequence may be missing.
        print("Warning: Missing data where FS expected or missing start. Internal stop may not be captured.")
        df_nd.loc[count, "ID"] = seq.id
    else:
#Confirm that methionine is present as first codon.
        if seqStringNH[1:3] == "TG":
            #Split sequences into codons and look for the first stop.
            codons = [(seqStringNH[i:i+3]) for i in range(0, len(seqStringNH), 3)]
            inFrameStopCodons = [i for i in codons if any(j in i for j in stops)]
            stopCodonIndeces = [i for i, x in enumerate(codons) if any(thing in x for thing in stops)]
         
#In case some seqs are missing the C-terminal region.
            if len(stopCodonIndeces) != 0:      
                firstStopLoc = stopCodonIndeces[0] * 3
            
#---BEGIN ANALYSES FOR IN-FRAME STOP SEQUENCES.---#        
            if len(stopCodonIndeces) != 0 and stopCodonIndeces[0]+1 < len(codons):
                seqStringNH1 = seqStringNH[firstStopLoc+1:]
                codons1 = [(seqStringNH1[i:i+3]) for i in range(0, len(seqStringNH1), 3)]
                inFrameStopCodons1 = [i for i in codons1 if any(j in i for j in stops)]
                stopCodonIndeces1 = [i for i, x in enumerate(codons1) if any(thing in x for thing in stops)]
                
                if len(stopCodonIndeces1) != 0:   
                    firstStopLoc1 = stopCodonIndeces1[0] * 3
                    
#Check to see if a second +1 frameshift is necessary to get the final protein.
                if len(stopCodonIndeces1) > 1:
                    part1 = seqStringToTranslate[0:firstStopLoc].translate(table)
                    part2 = seqStringToTranslate[firstStopLoc+1:firstStopLoc1].translate(table)
                    part3 = seqStringToTranslate[firstStopLoc1+1:].translate(table)
                    
                    fullTrans = "M" + str(part1[1:]) + str(part2) + str(part3)
                    
                    seqStringNH2 = seqStringNH1[firstStopLoc+1:]
                    codons2 = [(seqStringNH2[i:i+3]) for i in range(0, len(seqStringNH2), 3)]
                    inFrameStopCodons2 = [i for i in codons2 if any(j in i for j in stops)]
                    stopCodonIndeces2 = [i for i, x in enumerate(codons2) if any(thing in x for thing in stops)]
                    if len(stopCodonIndeces2) != 0:   
                        firstStopLoc2 = stopCodonIndeces2[0] * 3
#Give up trying to translate the sequence if it can't be translated by two +1 frameshifts.
                    if len(stopCodonIndeces2) > 1:
                        continue
#This is for strains with just one +1 frameshift
                else:
                    peptide = seqStringToTranslate[0:firstStopLoc].translate(table)
                    mainProt = seqStringToTranslate[firstStopLoc+1:].translate(table)
                    fullTrans = "M" + str(peptide[1:]) + str(mainProt)
                
                if "GGQ" in fullTrans:
                    #print(" ")
                    #print(seqStringNH)
                    df.loc[count, "in_frame_stop?"] = "yes"
                    df.loc[count, "internal_stop_index"] = firstStopLoc
                    df.loc[count, "SD_seq"] = seqStringNH[firstStopLoc-12:firstStopLoc-6]
                    df.loc[count, "FS_seq"] = seqStringNH[firstStopLoc-3:firstStopLoc]
                    df.loc[count, "internal_stop"] = seqStringNH[firstStopLoc:firstStopLoc+3]
                    df.loc[count, "nuc_after_stop"] = seqStringNH[firstStopLoc+3]
                    df.loc[count, "FS_context_seq"] = seqStringNH[firstStopLoc-18:firstStopLoc+11]
                    df.loc[count, "prfB_AA_seq"] = str(fullTrans)
#If there's not a GGQ motif in the translated sequence, it's likely inaccurate. Bin it.
                else:
                    continue
            
#---BEGIN ANALYSES FOR NONSTOP SEQUENCES.---#        
            else:
                #print("No stop codon.")
                transORF = seqStringToTranslate.translate(table)
                
                if "GGQ" in fullTrans:
                    df.loc[count, "in_frame_stop?"] = "no"
                    df.loc[count, "internal_stop_index"] = "NA"
                    df.loc[count, "SD_seq"] = seqString[expectedStop-12:expectedStop-6]
                    df.loc[count, "FS_seq"] = seqString[expectedStop-3:expectedStop]
                    df.loc[count, "internal_stop"] = seqString[expectedStop:expectedStop+3]
                    df.loc[count, "nuc_after_stop"] = seqString[expectedStop+3]
                    df.loc[count, "FS_context_seq"] = seqStringNH[expectedStop-18:expectedStop+11]
                    df.loc[count, "prfB_AA_seq"] = "M" + str(transORF[1:])
                    df.loc[count, "notes"] = "All extracted sequences come from expected locations based on the alignment file."

#Input general info about sequence.
            df.loc[count, "ID"] = seq.id
            df.loc[count, "seq_description"] = seq.description
            df.loc[count, "start"] = seqStringNH[0:3]
            df.loc[count, "prfB_DNA_seq"] = seqStringNH
            df.loc[count, "prfB_DNA_seq_alignment"] = seqString
            df.loc[count, "terminal_stop"] = seqStringNH[len(seqStringNH)-3:]

    count += 1

df_nd = df_nd.loc[(df_nd!=0).any(axis=1)]
df = df[df["in_frame_stop?"] != 0]

df.to_csv(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\FS_data_8_5_24.csv")
