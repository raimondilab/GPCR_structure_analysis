#!/usr/bin/python3

'''This module generates 3 figures for the interface analysis. One using all the structures, one using all
representative structures and one averaging all the structures representing the same GPCR-Galpha pair'''

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
import seaborn as sns
sns.set(color_codes=True, style="ticks", rc={'axes.facecolor':'white'}, font_scale=1.2)
colors = {"Gs": 'r', "Gi/o": 'b'}

def draw_violinplot(df, ax, data):

    '''Draw a violinplot and swarmplor with a p-value bar'''

    sns.violinplot(ax=ax, x=df["class"], y=df[data], palette=colors)
    sns.swarmplot(ax=ax, x=df["class"], y=df[data], color="black")
    pval = sp.mannwhitneyu(df[df["class"] == "Gi/o"][data], df[df["class"] == "Gs"][data])[1]
    x1, x2 = 0, 1
    y, ymin = max(np.array(df[data])), min(np.array(df[data]))
    step = (y-ymin)/5 # Characteristic length used to rescale the plot
    y += step
    ymin -= step
    h = step/5
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='gray')
    ax.text(0.5, y+1.5*h, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
    ax.set_ylim((ymin, y+step*0.8))

def draw_figure(df, filename, energy="Interface_energy", ratio="dG/dSASA", dSASA="dSASA"):

    '''Draw a figure with the violinplot of interface energy, dSASA and dG/dSASA'''

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    draw_violinplot(df, ax[0], energy)
    ax[0].set_ylabel("Interface energy (REU)")
    draw_violinplot(df, ax[1], dSASA)
    ax[1].set_ylabel("$\Delta$SASA (Å$^2$)")
    draw_violinplot(df, ax[2], ratio)
    ax[2].set_ylabel("$\Delta$SASA/$\Delta$G")
    fig.tight_layout()
    plt.savefig("../GPCR_experimental_structures/figures/"+filename)
    plt.clf()

total = pd.read_csv("GPCR_structs_clean.tsv", sep='\t')

total['class'] = total['Gprotein_Gene_name'].astype(str).str[:4]  # Divide Galphas by class
total = total[total["class"].isin(("GNAS", "GNAI", "GNAO", "GNAT"))]  # Delete underrepresented classes
total = total[total["Class"] == "classA"]
mask = total['Gprotein_Gene_name'] == "GNAS"
total.loc[mask, 'class'] = "Gs"
mask = np.logical_not(mask)  # Merge GNAI and GNAO classes
total.loc[mask, 'class'] = "Gi/o"
print("Number of Gi/o structures: ", len(total[total["class"] == "Gi/o"]))
print("Number of Gs structures: ", len(total[total["class"] == "Gs"]))
draw_figure(total, "binding_energy_total_cA.svg")

# Select representative structures: choose the one with the best resolution,
# in case of parity the one with the highest sequence coverage.
representative = total.sort_values(by=["Resolution", "Gprotein_residues", "Receptor_residues"], \
     key=lambda x: 10**9*total["Resolution"]-10**4*total["Gprotein_residues"]-total["Receptor_residues"])
representative = representative.drop_duplicates(subset=['Receptor_Gene_name', 'Gprotein_Gene_name'])
print("Number of Gi/o complexes: ", len(representative[representative["class"] == "Gi/o"]))
print("Number of Gs complexes: ", len(representative[representative["class"] == "Gs"]))
draw_figure(representative, "binding_energy_representative_cA.svg")

# Find the mean binding energy of all the structures of the same complex
condensed = total
condensed['mean_energy'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['Interface_energy'].transform('mean')
condensed['mean_dG/dSASA'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['dG/dSASA'].transform('mean')
condensed['mean_dSASA'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['dSASA'].transform('mean')
condensed = condensed.drop_duplicates(subset=['Receptor_Gene_name', 'Gprotein_Gene_name'])
draw_figure(condensed, "binding_energy_mean_cA.svg", energy="mean_energy", ratio='mean_dG/dSASA', dSASA='mean_dSASA')
