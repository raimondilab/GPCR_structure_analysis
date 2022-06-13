#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from functools import reduce
import seaborn as sns

fil=sys.argv[1]

df=pd.read_csv("./network/files/Links_struct_"+fil+".tsv", comment="#", sep="\t")
df1=pd.read_csv("./network/files/Nodes_struct_"+fil+".tsv", comment="#", sep="\t")

#degree
gs = df[df.cont_gs != 0]
gio = df[df.cont_gio != 0]
gq11 = df[df.cont_gq11 != 0]
#gpcr
degree_gpcr_gs=gs['SS1'].value_counts().reset_index()
degree_gpcr_gio=gio['SS1'].value_counts().reset_index()
degree_gpcr_gq11=gq11['SS1'].value_counts().reset_index()

data_frames = [degree_gpcr_gs, degree_gpcr_gio,degree_gpcr_gq11]

degree_gpcr = reduce(lambda  left,right: pd.merge(left,right,on=['index'],
                                            how='outer'), data_frames).fillna(0)


degree_gpcr.columns = ['GPCR', 'Gs','Gio','Gq11']
degree_gpcr.to_csv("./network/files/degree_gpcr_"+fil+".tsv",index=None,header=True,sep='\t')

#gprotein
degree_gprot_gs=gs['SS2'].value_counts().reset_index()
degree_gprot_gio=gio['SS2'].value_counts().reset_index()
degree_gprot_gq11=gq11['SS2'].value_counts().reset_index()
data_frames1 = [degree_gprot_gs, degree_gprot_gio,degree_gprot_gq11]
degree_gprot= reduce(lambda  left,right: pd.merge(left,right,on=['index'],
                                            how='outer'), data_frames1).fillna(0)

degree_gprot.columns = ['Gprotein', 'Gs','Gio','Gq11']
degree_gprot.to_csv("./network/files/degree_gprot_"+fil+".tsv",index=None,header=True,sep='\t')

#gpcr plot

degree_gpcr.set_index(['GPCR'],inplace=True)
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
fig, axes = plt.subplots(figsize=(15,15), facecolor=facecolor, ncols=2, sharey=True)
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
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=dgp,
    x=0, y="GPCR", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Node degree",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.02,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/degree_gpcr_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot degree
degree_gprot.set_index(['Gprotein'],inplace=True)
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(15,27), facecolor=facecolor, ncols=2, sharey=True)
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
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=dgp1,
    x=0, y="Gprotein", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("Gprotein structure",fontsize=50)
plt.ylabel("Node degree",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.1,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/degree_gprot_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')





###################################################################
###Node radius 
radius_gpcr= df.merge(df1, left_on=['SS1'],right_on=['struct'],how='left')
radius_gpcr=radius_gpcr.drop(radius_gpcr.columns[[0,1,2,3,4,5,6,7]], axis=1).drop_duplicates()
radius_gpcr.set_index(['struct'],inplace=True)
radius_gpcr=((radius_gpcr/radius_gpcr.values.max())*100).round(0)
radius_gpcr.columns = ['Gs','Gio','Gq11']
radius_gpcr.to_csv('./network/files/radius_gpcr_'+fil+'.tsv',index=True,header=True,sep='\t')

#GPROTEIN radius
radius_gprot= df.merge(df1, left_on=['SS2'],right_on=['struct'],how='left')
radius_gprot=radius_gprot.drop(radius_gprot.columns[[0,1,2,3,4,5,6,7]], axis=1).drop_duplicates()
radius_gprot.set_index(['struct'],inplace=True)
radius_gprot=((radius_gprot/radius_gprot.values.max())*100).round(0)
radius_gprot.columns = ['Gs','Gio','Gq11']
radius_gprot.to_csv('./network/files/radius_gprot_'+fil+'.tsv',index=True,header=True,sep='\t')


#gpcr plot
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(20,15), facecolor=facecolor, ncols=2, sharey=True)
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
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=rgp,
    x=0, y="struct", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Node radius",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.xticks(np.arange(0,110,10))

plt.legend(bbox_to_anchor=(0.02,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/radius_gpcr_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




#gprot plot
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(20,24), facecolor=facecolor, ncols=2, sharey=True)
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
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=rgp1,
    x=0, y="struct", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("Gprotein structure",fontsize=50)
plt.xlabel("Node radius",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(0.1,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.xticks(np.arange(0,110,10))

plt.savefig('./network/plots/radius_gprot_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

##########################################
###overall edge 'conservation
gs1 = df[df.fraction_gs != 0]
gio1 = df[df.fraction_gio != 0]
gq111 = df[df.fraction_gq11 != 0]
#gpcr
edge_gpcr_gs=gs1.groupby(['SS1'])['fraction_gs'].mean().reset_index()
edge_gpcr_gio=gio1.groupby(['SS1'])['fraction_gio'].mean().reset_index()
edge_gpcr_gq11=gq111.groupby(['SS1'])['fraction_gio'].mean().reset_index()

data_frame = [edge_gpcr_gs, edge_gpcr_gio,edge_gpcr_gq11]

edge_gpcr = reduce(lambda  left,right: pd.merge(left,right,on=['SS1'],
                                            how='outer'), data_frame).fillna(0)

edge_gpcr.columns = ['GPCR', 'Gs','Gio','Gq11']
edge_gpcr.set_index(['GPCR'],inplace=True)
edge_gpcr=(edge_gpcr*100).round(0)
edge_gpcr.to_csv('./network/files/edge_gpcr_'+fil+'.tsv',index=True,header=True,sep='\t')

#gprotein
edge_gprot_gs=gs1.groupby(['SS2'])['fraction_gs'].mean().reset_index()
edge_gprot_gio=gio1.groupby(['SS2'])['fraction_gio'].mean().reset_index()
edge_gprot_gq11=gq111.groupby(['SS2'])['fraction_gio'].mean().reset_index()

data_fram = [edge_gprot_gs, edge_gprot_gio,edge_gprot_gq11]

edge_gprot = reduce(lambda  left,right: pd.merge(left,right,on=['SS2'],
                                            how='outer'), data_fram).fillna(0)
edge_gprot.columns = ['Gprotein', 'Gs','Gio','Gq11']

edge_gprot.set_index(['Gprotein'],inplace=True)

edge_gprot=(edge_gprot*100).round(0)

edge_gprot.to_csv('./network/files/edge_gprot_'+fil+'.tsv',index=True,header=True,sep='\t')

#gpcr plot
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(22,15), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(edge_gpcr.index, edge_gpcr['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(edge_gpcr.index, edge_gpcr['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=edge_gpcr.index, yticklabels=edge_gpcr.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Edge Conservation(%)',fontsize=55)
axes[0].set_ylabel('GPCR position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0,edge_gpcr['Gio'].max()+5,5))
axes[1].xaxis.set_ticks(np.arange(0,edge_gpcr['Gio'].max()+5,5))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/edge_gpcr_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

#edge_gpcr.reset_index(['GPCR'],inplace=True)
egp=edge_gpcr.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,15))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=egp,
    x=0, y="GPCR", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("GPCR structure",fontsize=50)
plt.xlabel("Edge conservation(%)",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.xticks(np.arange(0,80,10))

plt.legend(bbox_to_anchor=(0.02,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/edge_gpcr_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')



#gprot degree
facecolor = '#eaeaf2'
color_red = '#c24a4d'
color_blue = '#4c72b0'
fig, axes = plt.subplots(figsize=(19,24), facecolor=facecolor, ncols=2, sharey=True)
fig.tight_layout()
axes[0].barh(edge_gprot.index, edge_gprot['Gs'], align='center', color=color_red, zorder=10)
axes[0].set_title('Gs', fontsize=65, pad=10, color=color_red)
axes[1].barh(edge_gprot.index, edge_gprot['Gio'], align='center', color=color_blue, zorder=10)
axes[1].set_title('Gio', fontsize=65, pad=10, color=color_blue)
axes[0].invert_xaxis() 
plt.gca().invert_yaxis()
axes[0].set(yticks=edge_gprot.index, yticklabels=edge_gprot.index)
axes[0].yaxis.tick_left()
axes[0].tick_params(axis='y', which='both', labelsize=55)
axes[0].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].tick_params(axis='x', which='both', labelsize=55,rotation=90)
axes[1].set_xlabel('Edge Conservation(%)',fontsize=55)

axes[0].set_ylabel('Gprotein position',fontsize=55)
axes[0].xaxis.set_ticks(np.arange(0, edge_gprot['Gs'].max()+5,5))
axes[1].xaxis.set_ticks(np.arange(0,edge_gprot['Gs'].max()+5 ,5))
axes[1].tick_params(axis='y', colors='#c24a4d')
plt.subplots_adjust(wspace=0, top=0.85, bottom=0.1, left=0.18, right=0.95)
axes[1].xaxis.set_label_coords(0, -0.1)
plt.savefig('./network/plots/edge_gprot_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')



egpr=edge_gprot.stack().reset_index()
sns.set(style='ticks')
plt.figure(figsize=(13,22))
palette ={"Gs": '#c24a4d', "Gio": '#4c72b0', "Gq11": '#008000'}
g = sns.barplot(
    data=egpr,
    x=0, y="Gprotein", hue="level_1"
   , palette=palette, alpha=1)
plt.ylabel("Gprotein structure",fontsize=50)
plt.xlabel("Edge conservation(%)",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.xticks(np.arange(0,50,10))
plt.legend(bbox_to_anchor=(0.02,1),fontsize=40,title='GPCR Class',title_fontsize=40,ncol=3,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./network/plots/edge_gprot_full_' +fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


