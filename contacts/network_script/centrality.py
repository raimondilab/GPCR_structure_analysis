#!/data/SW/anaconda3/envs/myenv/bin/python


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from functools import reduce
import seaborn as sns


df=pd.read_csv("all_gpcr.csv", comment="#", sep=",")
df1=pd.read_csv("all_gprot.csv", comment="#", sep=",")


#gpcr plot
sns.set(style='white')
df.set_index(['Nodes'],inplace=True)
facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
fig, axes = plt.subplots(figsize=(15,15), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(df.index, df['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(df.index, df['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=df.index, yticklabels=df.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].set_xlabel('Centrality betweenness',fontsize=55)
axes[0].set_ylabel('GPCR structure',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, 1.1,0.1))
axes[1].xaxis.set_ticks(np.arange(0,1.1 ,0.1))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.13)
plt.savefig('cent_gpcr.svg',bbox_inches='tight', dpi=600, format='svg')

dgp=df.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=dgp,
    x=0, y="Nodes", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Centrality betweenness",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
g.xaxis.set_ticks(np.arange(0,1.1 ,0.1))
plt.legend(bbox_to_anchor=(0.02,1),fontsize=40,title='Gprotein family',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('cent_gpcr_full.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot degree
sns.set(style='ticks')

df1.set_index(['Nodes'],inplace=True)
facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(15,27), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(df1.index, df1['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(df1.index, df1['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=df1.index, yticklabels=df1.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].set_xlabel('Centrality betweenness',fontsize=55)
axes[0].set_ylabel('Gprotein structure',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, 1.1,0.1))
axes[1].xaxis.set_ticks(np.arange(0,1.1 ,0.1))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('cent_gprot.svg',bbox_inches='tight', dpi=600, format='svg')


dgp1=df1.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=dgp1,
    x=0, y="Nodes", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Centrality betweenness",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
g.xaxis.set_ticks(np.arange(0,1.1 ,0.1))
plt.legend(bbox_to_anchor=(1,1.16),fontsize=40,title='Gprotein family',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('cent_gprot_full.svg',bbox_inches='tight', dpi=600, format='svg')



#####classA
df=pd.read_csv("all_gpcr_classA.csv", comment="#", sep=",")
df1=pd.read_csv("all_gprot_classA.csv", comment="#", sep=",")


#gpcr plot
sns.set(style='ticks')
df.set_index(['Nodes'],inplace=True)
facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
fig, axes = plt.subplots(figsize=(15,15), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(df.index, df['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(df.index, df['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=df.index, yticklabels=df.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].set_xlabel('Centrality betweenness',fontsize=55)
axes[0].set_ylabel('GPCR structure',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, 1.1,0.1))
axes[1].xaxis.set_ticks(np.arange(0,1.1 ,0.1))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.13)
plt.savefig('cent_gpcr_classA.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot degree
sns.set(style='ticks')

df1.set_index(['Nodes'],inplace=True)
facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(15,27), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(df1.index, df1['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(df1.index, df1['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=df1.index, yticklabels=df1.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=50,rotation=90)
axes[1].set_xlabel('Centrality betweenness',fontsize=55)
axes[0].set_ylabel('Gprotein structure',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, 1.1,0.1))
axes[1].xaxis.set_ticks(np.arange(0,1.1 ,0.1))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('cent_gprot_classA.svg',bbox_inches='tight', dpi=600, format='svg')

