#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
import seaborn as sns; sns.set(color_codes=True)

df = pd.read_csv('../GPCR_experimental_structures/binding_energy.tsv', sep='\t')
df = df[df.Interface_energy<0] # Delete unrelaiable positive binding energies
df['class']=df['Gprotein_Gene_name'].astype(str).str[:4] # Divide Galphas by class
df=df[df["class"].isin(("GNAS","GNAI","GNAO"))] # Delete underrepresented classes
mask=df['Gprotein_Gene_name']=="GNAS"
df.loc[mask, 'class'] = "Gs"
mask=(np.logical_or(df['Gprotein_Gene_name'].astype(str).str[:4]=="GNAI",df['Gprotein_Gene_name'].astype(str).str[:4]=="GNAO")) # Merge GNAI and GNAO classes
df.loc[mask, 'class'] = "Gi/o"
print("Number of Gi/o structures: ",len(df[df["class"]=="Gi/o"]))
print("Number of Gs structures: ",len(df[df["class"]=="Gs"]))

# First violin plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})

ax = sns.violinplot(x=df["class"],y=df["Interface_energy"])
pval=sp.ttest_ind(df[df["class"]=="Gi/o"].Interface_energy,df[df["class"]=="Gs"].Interface_energy)[1]
x1, x2 = 0, 1
y, ymin,h, col = max(np.array(df["Interface_energy"])) + 15,min(np.array(df["Interface_energy"])) -20, 4, 'gray'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text(0.5, y+h+2, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
plt.ylim((ymin,y+18))
plt.savefig('../GPCR_experimental_structures/figures/binding_energy_total.svg')
plt.clf()

# Select representative structures (choose the one with the best resolution, in case of parity the one with the highest sequence coverage)
representative=df.sort_values(by=["Resolution","Gprotein_residues","Receptor_residues"],key=lambda x: 10**9*df["Resolution"]-df["Receptor_residues"]*df["Gprotein_residues"])
representative=representative.drop_duplicates(subset=['Receptor_Uniprot_AC', 'Gprotein_Uniprot_AC'])
print("Number of Gi/o complexes: ",len(representative[representative["class"]=="Gi/o"]))
print("Number of Gs complexes: ",len(representative[representative["class"]=="Gs"]))

# Second violin plot
ax = sns.violinplot(x=representative["class"],y=representative["Interface_energy"])
pval=sp.ttest_ind(representative[representative["class"]=="Gi/o"].Interface_energy,representative[representative["class"]=="Gs"].Interface_energy)[1]
x1, x2 = 0, 1
y, ymin,h, col = max(np.array(representative["Interface_energy"])) + 15,min(np.array(representative["Interface_energy"])) -20, 4, 'gray'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text(0.5, y+h+2, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
plt.ylim((ymin,y+18))
plt.savefig('../GPCR_experimental_structures/figures/binding_energy_representative.svg')
plt.clf()

# Find the mean binding energy of all the structures of the same complex
condensed=df
condensed['mean_energy'] = condensed.groupby(['Receptor_Uniprot_AC', 'Gprotein_Uniprot_AC'])['Interface_energy'].transform('mean')
condensed=condensed.drop_duplicates(subset=['Receptor_Uniprot_AC', 'Gprotein_Uniprot_AC'])

# Third violin plot
ax = sns.violinplot(x=condensed["class"],y=condensed["mean_energy"])
pval=sp.ttest_ind(condensed[condensed["class"]=="Gi/o"].mean_energy,condensed[condensed["class"]=="Gs"].mean_energy)[1]
x1, x2 = 0, 1
y, ymin,h, col = max(np.array(condensed["mean_energy"])) + 15,min(np.array(condensed["mean_energy"])) -20, 4, 'gray'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text(0.5, y+h+2, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
plt.ylim((ymin,y+18))
plt.savefig('../GPCR_experimental_structures/figures/binding_energy_mean.svg')
