import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from functools import reduce
import seaborn as sns
import networkx as nx
#################################### Calculate network statistics degree,radius, and centrality

fil=sys.argv[1]

df=pd.read_csv("./network/files/Links_struct_"+fil+".tsv", comment="#", sep="\t")
df1=pd.read_csv("./network/files/Nodes_struct_"+fil+".tsv", comment="#", sep="\t")

#degree
gs = df[df.cont_gs != 0]
gio = df[df.cont_gio != 0]
gq11 = df[df.cont_gq11 != 0]
g1213 = df[df.cont_g1213 != 0]

#gpcr
degree_gpcr_gs=gs['SS1'].value_counts().reset_index()
degree_gpcr_gio=gio['SS1'].value_counts().reset_index()
degree_gpcr_gq11=gq11['SS1'].value_counts().reset_index()
degree_gpcr_g1213=g1213['SS1'].value_counts().reset_index()

data_frames = [degree_gpcr_gs, degree_gpcr_gio,degree_gpcr_gq11,degree_gpcr_g1213]

degree_gpcr = reduce(lambda  left,right: pd.merge(left,right,on=['index'],
                                            how='outer'), data_frames).fillna(0)


degree_gpcr.columns = ['GPCR', 'Gs','Gio','Gq11','G1213']
degree_gpcr.to_csv("./statistics/network_stat/degree_gpcr_"+fil+".tsv",index=None,header=True,sep='\t')

#gprotein
degree_gprot_gs=gs['SS2'].value_counts().reset_index()
degree_gprot_gio=gio['SS2'].value_counts().reset_index()
degree_gprot_gq11=gq11['SS2'].value_counts().reset_index()
degree_gprot_g1213=g1213['SS2'].value_counts().reset_index()

data_frames1 = [degree_gprot_gs, degree_gprot_gio,degree_gprot_gq11,degree_gprot_g1213]
degree_gprot= reduce(lambda  left,right: pd.merge(left,right,on=['index'],
                                            how='outer'), data_frames1).fillna(0)

degree_gprot.columns = ['Gprotein', 'Gs','Gio','Gq11','G1213']

degree_gprot.to_csv("./statistics/network_stat/degree_gprot_"+fil+".tsv",index=None,header=True,sep='\t')

#gpcr plot

degree_gpcr.set_index(['GPCR'],inplace=True)
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
color_purple='#800080'

fig, axes = plt.subplots(figsize=(15,15), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(degree_gpcr.index, degree_gpcr['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(degree_gpcr.index, degree_gpcr['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=degree_gpcr.index, yticklabels=degree_gpcr.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Node Degree',fontsize=55)
axes[0].set_ylabel('GPCR position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, degree_gpcr.values.max(),2))
axes[1].xaxis.set_ticks(np.arange(0,degree_gpcr.values.max() ,2))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/degree_gpcr_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

dgp=degree_gpcr.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=dgp,
    x=0, y="GPCR", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Node degree",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.01,1.01),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=2,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/degree_gpcr_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot degree
degree_gprot.set_index(['Gprotein'],inplace=True)
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(15,27), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(degree_gprot.index, degree_gprot['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(degree_gprot.index, degree_gprot['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=degree_gprot.index, yticklabels=degree_gprot.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Node Degree',fontsize=55)
axes[0].set_ylabel('Gprotein position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, degree_gprot.values.max(),2))
axes[1].xaxis.set_ticks(np.arange(0,degree_gprot.values.max() ,2))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/degree_gprot_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


dgp1=degree_gprot.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(15,30))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=dgp1,
    x=0, y="Gprotein", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("Gprotein structure",fontsize=50)
plt.ylabel("Node degree",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.1,1),fontsize=40,title='Gfamily',title_fontsize=40,ncol=2,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/degree_gprot_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')





###################################################################
###Node radius 
radius_gpcr= df.merge(df1, left_on=['SS1'],right_on=['struct'],how='left')
radius_gpcr=radius_gpcr.drop(radius_gpcr.columns[[0,1,2,3,4,5,6,7,8,9]], axis=1).drop_duplicates()
radius_gpcr.set_index(['struct'],inplace=True)
radius_gpcr=((radius_gpcr/radius_gpcr.values.max())*100).round(0)
radius_gpcr.columns = ['Gs','Gio','Gq11','G1213']
radius_gpcr.to_csv('./statistics/network_stat/radius_gpcr_'+fil+'.tsv',index=True,header=True,sep='\t')

#GPROTEIN radius
radius_gprot= df.merge(df1, left_on=['SS2'],right_on=['struct'],how='left')
radius_gprot=radius_gprot.drop(radius_gprot.columns[[0,1,2,3,4,5,6,7,8,9]], axis=1).drop_duplicates()
radius_gprot.set_index(['struct'],inplace=True)
radius_gprot=((radius_gprot/radius_gprot.values.max())*100).round(0)
radius_gprot.columns = ['Gs','Gio','Gq11','G1213']
radius_gprot.to_csv('./statistics/network_stat/radius_gprot_'+fil+'.tsv',index=True,header=True,sep='\t')


#gpcr plot
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(20,15), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(radius_gpcr.index, radius_gpcr['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(radius_gpcr.index, radius_gpcr['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=radius_gpcr.index, yticklabels=radius_gpcr.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Node Radius',fontsize=55)
axes[0].set_ylabel('GPCR position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, radius_gpcr.values.max()+10,10))
axes[1].xaxis.set_ticks(np.arange(0,radius_gpcr.values.max()+10,10))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/radius_gpcr_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

rgp=radius_gpcr.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=rgp,
    x=0, y="struct", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Node radius",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.xticks(np.arange(0,110,10))

plt.legend(bbox_to_anchor=(0.02,1.01),fontsize=40,title='Gfamily',title_fontsize=40,ncol=2,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/radius_gpcr_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot plot
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(20,24), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(radius_gprot.index, radius_gprot['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(radius_gprot.index, radius_gprot['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=radius_gprot.index, yticklabels=radius_gprot.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Node Radius',fontsize=55)

axes[0].set_ylabel('Gprotein position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, radius_gprot.values.max()+10,10))
axes[1].xaxis.set_ticks(np.arange(0,radius_gprot.values.max()+10 ,10))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/radius_gprot_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')



rgp1=radius_gprot.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(15,30))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=rgp1,
    x=0, y="struct", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("Gprotein structure",fontsize=50)
plt.xlabel("Node radius",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.1,1.01),fontsize=40,title='Gfamily',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.xticks(np.arange(0,110,10))
plt.savefig('./network/plots/radius_gprot_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')



########################################centrality
#make graphs
s = nx.Graph()
s1 = nx.Graph()
s2 = nx.Graph()
s3 = nx.Graph()
s.add_nodes_from(degree_gpcr_gs['index'].tolist())
s1.add_nodes_from(degree_gpcr_gio['index'].tolist())
s2.add_nodes_from(degree_gpcr_gq11['index'].tolist())
s3.add_nodes_from(degree_gpcr_g1213['index'].tolist())
s.add_nodes_from(degree_gprot_gs['index'].tolist())
s1.add_nodes_from(degree_gprot_gio['index'].tolist())
s2.add_nodes_from(degree_gprot_gq11['index'].tolist())
s3.add_nodes_from(degree_gprot_g1213['index'].tolist())
s.add_edges_from([tuple(r) for r in gs[['SS1', 'SS2']].to_numpy()])
s1.add_edges_from([tuple(r) for r in gio[['SS1', 'SS2']].to_numpy()])
s2.add_edges_from([tuple(r) for r in gq11[['SS1', 'SS2']].to_numpy()])
s3.add_edges_from([tuple(r) for r in g1213[['SS1', 'SS2']].to_numpy()])
#calculate centrality
btwn =  nx.betweenness_centrality(s)
btwn1 =  nx.betweenness_centrality(s1)
btwn2 =  nx.betweenness_centrality(s2)
btwn3 =  nx.betweenness_centrality(s3)


bt=pd.DataFrame(btwn.items(), columns=['Nodes', 'Gs'])
bt1=pd.DataFrame(btwn1.items(), columns=['Nodes', 'Gio'])
bt2=pd.DataFrame(btwn2.items(), columns=['Nodes', 'Gq11'])
bt3=pd.DataFrame(btwn3.items(), columns=['Nodes', 'G1213'])
data_frames = [bt, bt1,bt2,bt3]

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['Nodes'],
                                            how='outer'), data_frames)
df_merged=df_merged.fillna(0)
df_merged.set_index(['Nodes'],inplace=True)
df_merged=df_merged.loc[~(df_merged==0).all(axis=1)]
df_merged=df_merged.round(3)
gpcr_cent=df_merged[df_merged.index.str.contains(r'TM|ICL|H8|Ct|Nt')]
gprot_cent=df_merged[~df_merged.isin(gpcr_cent)].dropna()
gpcr_cent.to_csv("./statistics/network_stat/cent_gpcr_"+fil+".tsv",index=None,header=True,sep='\t')
gprot_cent.to_csv("./statistics/network_stat/cent_gprot_"+fil+".tsv",index=None,header=True,sep='\t')

#gpcr plot
sns.set(style='white')
facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
fig, axes = plt.subplots(figsize=(15,15), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(gpcr_cent.index, gpcr_cent['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(gpcr_cent.index, gpcr_cent['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=gpcr_cent.index, yticklabels=gpcr_cent.index)
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
plt.savefig('./network/plots/cent_gpcr_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

dgp=gpcr_cent.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=dgp,
    x=0, y="Nodes", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Centrality betweenness",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
g.xaxis.set_ticks(np.arange(0,1.1 ,0.1))
plt.legend(bbox_to_anchor=(0.02,1.01),fontsize=40,title='Gprotein family',title_fontsize=40,ncol=2,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/cent_gpcr_full_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

#gprot 
sns.set(style='ticks')

facecolor = 'white'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(15,27), ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(gprot_cent.index, gprot_cent['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(gprot_cent.index, gprot_cent['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=gprot_cent.index, yticklabels=gprot_cent.index)
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
plt.savefig('./network/plots/cent_gprot_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


dgp1=gprot_cent.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000','G1213':'#800080'}
g = sns.barplot(
    data=dgp1,
    x=0, y="Nodes", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Centrality betweenness",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
g.xaxis.set_ticks(np.arange(0,1.1 ,0.1))
plt.legend(bbox_to_anchor=(0.95,1.23),fontsize=40,title='Gprotein family',title_fontsize=40,ncol=2,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/cent_gprot_full_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')
