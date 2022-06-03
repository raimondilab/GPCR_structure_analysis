#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sp
import numpy as np
import seaborn as sns; sns.set(color_codes=True)

def draw_violinplot(df,ax,data):
	sns.violinplot(ax=ax,x=df["class"],y=df[data])
	pval=sp.mannwhitneyu(df[df["class"]=="Gi/o"][data],df[df["class"]=="Gs"][data])[1]
	x1, x2 = 0, 1
	y, ymin = max(np.array(df[data])),min(np.array(df[data]))
	step=(y-ymin)/5
	y+=step
	ymin-=step
	h=step/5
	ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c='gray')
	ax.text(0.5, y+1.5*h, "P-val=%.1E" % pval, ha='center', va='bottom', color='black')
	ax.set_ylim((ymin,y+step*0.8))
	return
	
def draw_figure(df,filename,energy="Interface_energy",ratio="dG/dSASA",dSASA="dSASA"):
	fig, ax = plt.subplots(1, 3,figsize=(15,5))
	draw_violinplot(df,ax[0],energy)
	draw_violinplot(df,ax[1],dSASA)
	draw_violinplot(df,ax[2],ratio)
	fig.tight_layout()
	plt.savefig("../GPCR_experimental_structures/figures/"+filename)
	plt.clf()
	return

df = pd.read_csv("../GPCR_experimental_structures/binding_energy.tsv", sep='\t')

#df = df[df.Interface_energy<0] # Delete unrelaiable positive binding energies
df['class']=df['Gprotein_Gene_name'].astype(str).str[:4] # Divide Galphas by class
df=df[df["class"].isin(("GNAS","GNAI","GNAO"))] # Delete underrepresented classes
mask=df['Gprotein_Gene_name']=="GNAS"
df.loc[mask, 'class'] = "Gs"
mask=(np.logical_or(df['Gprotein_Gene_name'].astype(str).str[:4]=="GNAI",df['Gprotein_Gene_name'].astype(str).str[:4]=="GNAO")) # Merge GNAI and GNAO classes
df.loc[mask, 'class'] = "Gi/o"
#df["Interface_energy"]=df["Interface_energy"]*50
print("Number of Gi/o structures: ",len(df[df["class"]=="Gi/o"]))
print("Number of Gs structures: ",len(df[df["class"]=="Gs"]))

# First plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})

draw_figure(df,"binding_energy_total.svg")

# Select representative structures (choose the one with the best resolution, in case of parity the one with the highest sequence coverage)
representative=df.sort_values(by=["Resolution","Gprotein_residues","Receptor_residues"],key=lambda x: 10**9*df["Resolution"]-10**4*df["Receptor_residues"]-df["Gprotein_residues"])
representative=representative.drop_duplicates(subset=['Receptor_Gene_name', 'Gprotein_Gene_name'])
print("Number of Gi/o complexes: ",len(representative[representative["class"]=="Gi/o"]))
print("Number of Gs complexes: ",len(representative[representative["class"]=="Gs"]))
draw_figure(representative,"binding_energy_representative.svg")

# Find the mean binding energy of all the structures of the same complex
condensed=df
condensed['mean_energy'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['Interface_energy'].transform('mean')
condensed['mean_dG/dSASA'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['dG/dSASA'].transform('mean')
condensed['mean_dSASA'] = condensed.groupby(['Receptor_Gene_name', 'Gprotein_Gene_name'])['dSASA'].transform('mean')
condensed=condensed.drop_duplicates(subset=['Receptor_Gene_name', 'Gprotein_Gene_name'])
draw_figure(condensed,"binding_energy_mean.svg",energy="mean_energy",ratio='mean_dG/dSASA',dSASA='mean_dSASA')
