# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:19:12 2023

@author: cassp
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker

df = pd.read_csv(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Data\Bioinformatics\FS_data_8_5_24.csv")


###--- BEGIN MAKING PLOTS ---###

df_y = df[df["in_frame_stop?"] == "yes"]
df_y_FS = df_y["FS_context_seq"].dropna()
fsDf_y = df_y_FS.apply(lambda x: pd.Series(list(x)))
df_n = df[df["in_frame_stop?"] == "no"]
df_n_FS = df_n["FS_context_seq"].dropna()
fsDf_n = df_n_FS.apply(lambda x: pd.Series(list(x)))

totals_y = pd.DataFrame(0, index = range(28), columns = ["Position", "A", "T", "G", "C"])
totals_y_prop = pd.DataFrame(0, index = range(28), columns = ["A", "T", "G", "C"])
totals_n = pd.DataFrame(0, index = range(28), columns = ["Position", "A", "T", "G", "C"])
totals_n_prop = pd.DataFrame(0, index = range(28), columns = ["A", "T", "G", "C"])

position = 0
while position < fsDf_y.shape[1]:
    short = fsDf_y.iloc[:, position]
    totals_y.loc[position, "Position"] = position
    totals_y.loc[position, "A"] = np.count_nonzero(short == "A")
    totals_y.loc[position, "T"] = np.count_nonzero(short == "T")
    totals_y.loc[position, "G"] = np.count_nonzero(short == "G")
    totals_y.loc[position, "C"] = np.count_nonzero(short == "C")
    
    totals_y_prop.loc[position, "A"] = totals_y.loc[position, "A"]/sum(totals_y.iloc[position, 1:5])
    totals_y_prop.loc[position, "T"] = totals_y.loc[position, "T"]/sum(totals_y.iloc[position, 1:5])
    totals_y_prop.loc[position, "G"] = totals_y.loc[position, "G"]/sum(totals_y.iloc[position, 1:5])
    totals_y_prop.loc[position, "C"] = totals_y.loc[position, "C"]/sum(totals_y.iloc[position, 1:5])
    
    position += 1

position = 0
while position < fsDf_n.shape[1]:
    short = fsDf_n.iloc[:, position]
    totals_n.loc[position, "Position"] = position
    totals_n.loc[position, "A"] = np.count_nonzero(short == "A")
    totals_n.loc[position, "T"] = np.count_nonzero(short == "T")
    totals_n.loc[position, "G"] = np.count_nonzero(short == "G")
    totals_n.loc[position, "C"] = np.count_nonzero(short == "C")
    
    totals_n_prop.loc[position, "A"] = totals_n.loc[position, "A"]/sum(totals_n.iloc[position, 1:5])
    totals_n_prop.loc[position, "T"] = totals_n.loc[position, "T"]/sum(totals_n.iloc[position, 1:5])
    totals_n_prop.loc[position, "G"] = totals_n.loc[position, "G"]/sum(totals_n.iloc[position, 1:5])
    totals_n_prop.loc[position, "C"] = totals_n.loc[position, "C"]/sum(totals_n.iloc[position, 1:5])
    
    position += 1


logo = logomaker.Logo(totals_y_prop, color_scheme = "classic")
logo.ax.set_ylabel('Proportion', fontsize=16)
logo.ax.set_xlabel('Position', fontsize=16)
logo.ax.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
plt.savefig(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Figures\FS_logo_8_5_24.png", dpi=600)

logo = logomaker.Logo(totals_n_prop, color_scheme = "classic")
logo.ax.set_ylabel('Proportion', fontsize=16)
logo.ax.set_xlabel('Position', fontsize=16)
logo.ax.tick_params(axis='both', which='major', labelsize=12)
plt.tight_layout()
plt.savefig(r"C:\Users\cassp\Box Sync\Feaga Lab\Cassidy Prince\prfB\Figures\no_FS_logo_8_5_24.png", dpi=600)
