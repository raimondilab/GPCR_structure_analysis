#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.patches as mpatches


df=pd.read_csv("cluster_ratio_rmsd.txt", comment="#", sep="\t")




color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
color_purple='#800080'
colors = [color_red,color_blue,color_green,color_purple,color_red,color_blue,color_green,color_purple,color_red,color_blue,color_green,color_purple]

#cluster1
sns.set(style='ticks')
plt.figure(figsize=(15,15))

g = sns.barplot(
    data=df,
    x="Plot", y="val", alpha=1,palette=colors,hue='Cluster', edgecolor = 'black')

plt.xlabel("Gprotein Family",fontsize=60)
plt.ylabel("Fraction of coupling",fontsize=60)
plt.xticks(fontsize=60,rotation=90)
plt.yticks(np.arange(0,1.1,0.1),fontsize=60)
plt.title('Clusters',fontsize=60)
hatches = ['\\','\\','\\','\\','+','+','+','+','-','-','-','-']
# Loop over the bars
for i,thisbar in enumerate(g.patches):
    # Set a different hatch for each bar
    thisbar.set_color(colors[i])
for i,thisbar in enumerate(g.patches):
    thisbar.set_edgecolor("black")
    thisbar.set_linewidth(3)
# Loop over the bars
for i,s in enumerate(g.patches):
    # Set a different hatch for each bar
    s.set_hatch(hatches[i])
clu1 = mpatches.Patch(edgecolor='black', facecolor='white',
                     hatch='\\', label='Cluster 1')
clu2 = mpatches.Patch(edgecolor='black', facecolor='white', label='Cluster 2',hatch='+')
clu3 = mpatches.Patch(edgecolor='black', facecolor='white', label='Cluster 3',hatch='-')

plt.legend(handles=[clu1,clu2,clu3],bbox_to_anchor=(0.66,0.98),fontsize=30,title='Gprotein Family',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0.1 ,framealpha=1)


plt.savefig('clusters_all_rmsd.pdf',bbox_inches='tight', dpi=600, format='pdf')

