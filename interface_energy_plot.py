#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
import seaborn as sns; sns.set(color_codes=True)

df = pd.read_csv('../GPCR_experimental_structures/binding_energy.tsv', sep='\t')
condensed=df[['Receptor_Uniprot_AC','Receptor_Gene_name','Receptor_Uniprot_ID','Gprotein_Uniprot_AC','Gprotein_Gene_name','Gprotein_Uniprot_ID','Interface_energy']]
condensed=condensed[condensed.Interface_energy<0] # Delete unrelaiable positive bing energies
condensed['class']=condensed['Gprotein_Gene_name'].astype(str).str[:4] # Divide Galphas by class
condensed=condensed[condensed["class"]!="GNAT"] # Delete GNAT
mask=condensed['Gprotein_Gene_name']=="GNAS"
condensed.loc[mask, 'class'] = "Gs"
mask=(np.logical_or(condensed['Gprotein_Gene_name'].astype(str).str[:4]=="GNAI",condensed['Gprotein_Gene_name'].astype(str).str[:4]=="GNAO"))
condensed.loc[mask, 'class'] = "Gi/o" # Merge GNAI and GNAO classes

condensed['mean_energy'] = condensed.groupby(['Receptor_Uniprot_AC', 'Gprotein_Uniprot_AC'])['Interface_energy'].transform('mean') # Average the binding energies of the same complexes
condensed=condensed.drop_duplicates(subset=['Receptor_Uniprot_AC', 'Gprotein_Uniprot_AC'])

# Violin plot

sns.set(style="ticks",rc={'axes.facecolor':'white'})

ax = sns.violinplot(x=condensed["class"],y=condensed["mean_energy"])
pval=sp.ttest_ind(condensed[condensed["class"]=="Gi/o"].mean_energy,condensed[condensed["class"]=="Gs"].mean_energy)[1]
x1, x2 = 0, 1
y, ymin,h, col = max(np.array(condensed["mean_energy"])) + 15,min(np.array(condensed["mean_energy"])) -20, 4, 'gray'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text(0.5, y+h+1.5, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
plt.ylim((ymin,0))
plt.savefig('../GPCR_experimental_structures/figures/binding_energy_mean.svg')
