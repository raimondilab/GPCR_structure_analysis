#!/data/SW/anaconda3/envs/myenv/bin/python

from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from statannotations.Annotator import Annotator
import sys 
import numpy as np
fs='\t'
fil=sys.argv[1]
#Gs vs Gi pairs
path='./plots/'
df=pd.read_csv(fil+".txt",header=None, comment="#", sep="\t")
df.columns=['PDB_ID','PDB_ID1','rmsd']

df = df.replace(r'\n',' ', regex=True) 
df['PDB_ID'] = df['PDB_ID'].astype(str).str.strip()
df['PDB_ID1'] = df['PDB_ID1'].astype(str).str.strip()

df1=pd.read_csv("GPCR_structs_clean_consensus.tsv", comment="#", sep="\t")
ds=df1[['PDB_ID','Receptor_Gene_name', 'Gprotein_Gene_name']].drop_duplicates()
ds['name']=ds['Receptor_Gene_name']+"_"+ ds['Gprotein_Gene_name']+"_"+ds['PDB_ID']

ds=ds[['PDB_ID','name']].drop_duplicates()
ds['PDB_ID'] = ds['PDB_ID'].astype(str)
ds['name'] = ds['name'].astype(str)
ds=ds.reset_index(drop=True)
mydict = dict(zip(ds.PDB_ID, ds.name))
df2=df.replace({"PDB_ID": mydict})
df2=df2.replace({"PDB_ID1": mydict})
df2[['GPCR1', 'Gprotein1', 'pdb1']] = df2['PDB_ID'].str.split('_', expand=True)
df2[['GPCR2', 'Gprotein2', 'pdb2']] = df2['PDB_ID1'].str.split('_', expand=True)

equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df2["gfam1"] = df2["Gprotein1"].map(equiv)
df2["gfam2"] = df2["Gprotein2"].map(equiv)
df2=df2.dropna()
df2.to_csv(fil+".tsv",index=None,header=True,sep='\t')

#draw
ds=df2[['PDB_ID', 'gfam1']].drop_duplicates()
ds1=df2[['PDB_ID1', 'gfam2']].drop_duplicates()
ds.columns=['PDB_ID','gfam']
ds1.columns=['PDB_ID','gfam']
ds2=pd.concat([ds,ds1]).drop_duplicates()
ds2.set_index('PDB_ID',inplace=True)
ds2 = ds2.pop("gfam")
ds2=ds2.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}
rs = ds2.map(lut)
#class
dc=df2[['PDB_ID', 'GPCR1']].drop_duplicates()
dc1=df2[['PDB_ID1', 'GPCR2']].drop_duplicates()
dc.columns=['PDB_ID','GPCR']
dc1.columns=['PDB_ID','GPCR']
dc2=pd.concat([dc,dc1]).drop_duplicates()
clas=pd.read_csv("./use_file/classification.txt", comment="#", sep="\t")
clas1 = dc2.merge(clas,left_on=['GPCR',],right_on=['GPCR'],how="left")
clas1=clas1[['PDB_ID', 'class']].drop_duplicates()
clas1.set_index('PDB_ID',inplace=True)
clas1 = clas1.pop("class")
clas1=clas1.astype('object')
clas1=clas1.astype('object')
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}
rs1 = clas1.map(lut1)

#couplings

gp=pd.read_csv("./use_file/coup_assay.txt", sep="\t",index_col=0)
gp1 = dc2.merge(gp,left_on=['GPCR'],right_on=['GPCR'],how="left")
gp1.set_index('PDB_ID',inplace=True)
gp1=gp1.drop(columns=['GPCR'])
rs3 = pd.concat([rs,rs1],axis=1)
rs3 = pd.concat([rs3,gp1],axis=1)
lut.update(lut1)
coup={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}
lut.update(coup) 
ab=pd.read_csv("./use_file/antibody.txt", sep="\t",index_col=0)
ab1={'Nb35 antibody':'#da70d6','Scfv16 antibody':'#b8860b','Other antibodies':'#808000'}
lut.update(ab1) 
rs3 = pd.concat([rs3,ab],axis=1)


#classA
clas1=clas1.reset_index()
mydict1 = dict(zip(clas1['PDB_ID'],clas1['class']))
df2["class1"] = df2["PDB_ID"].map(mydict1)
df2["class2"] = df2["PDB_ID1"].map(mydict1)
df3=df2.loc[((df2['class1'] =='classA') & (df2['class2'] =='classA'))] 

dh1=df2[['PDB_ID','PDB_ID1', 'rmsd']].drop_duplicates()
dh2=df2[['PDB_ID','PDB_ID1', 'rmsd']].drop_duplicates()
dh2.columns=['PDB_ID1','PDB_ID', 'rmsd']
dh3 = pd.concat([dh1,dh2]).drop_duplicates()

table1 = pd.pivot_table(dh3,values='rmsd', index=['PDB_ID'],columns=['PDB_ID1']).fillna(0)
da1=df3[['PDB_ID','PDB_ID1', 'rmsd']].drop_duplicates()
da2=df3[['PDB_ID','PDB_ID1', 'rmsd']].drop_duplicates()
da2.columns=['PDB_ID1','PDB_ID', 'rmsd']
da3 = pd.concat([da1,da2]).drop_duplicates()

table2 = pd.pivot_table(da3,values='rmsd', index=['PDB_ID'],columns=['PDB_ID1']).fillna(0)
mask=0
#all class
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})

g=sns.clustermap(table1,method='ward',mask=table1==mask,row_colors=rs3,col_colors=rs3,tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.04, 0.04),cmap='Reds',square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver',colors_ratio=(0.01,0.01),figsize=(27,27))
g.ax_heatmap.set_xlabel('Position', fontsize=25)
g.ax_heatmap.set_ylabel('GPCR', fontsize=25)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =16)
g.ax_heatmap.tick_params(axis='x', which='major', pad=2)
g.ax_heatmap.tick_params(axis='y', which='major', pad=2)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 16) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
#g.gs.update(top=0.1,right=0.05)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1.columns), linewidth=4, color='black')
g.ax_heatmap.hlines(y=len(table1.index), xmin=-0.5, xmax=len(table1.columns), linewidth=5, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1.index), linewidth=5, color='black')
g.ax_heatmap.vlines(x=len(table1.columns), ymin=0, ymax=len(table1.index), linewidth=5, color='black')
g.cax.set_position([1.1, 0.45, 0.02, 0.3]) 
#g.cax.set_yticklabels(g.cax.get_ymajorticklabels(), fontsize =22)
g.cax.tick_params(labelsize=22)
g.ax_row_colors.tick_params(labelsize=18)
g.ax_col_colors.tick_params(labelsize=18)

g.cax.set_title("RMSD",loc="center",pad=20,fontsize=22)
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.7,-0.1),fontsize=20,title_fontsize=20, bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=5,title='Gprotein family/GPCR class')


g.savefig(path+fil+'htmp_all.svg',bbox_inches='tight',dpi=300)

#classA
#all class
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})

g=sns.clustermap(table2,method='ward',mask=table2==mask,row_colors=rs3,col_colors=rs3,tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.04, 0.04),cmap='Reds',square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver',colors_ratio=(0.01,0.01),figsize=(17,17))
g.ax_heatmap.set_xlabel('Position', fontsize=45)
g.ax_heatmap.set_ylabel('GPCR', fontsize=45)
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =35)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 35) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.gs.update(top=3.5,right=3.5)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table2.columns), linewidth=5, color='black')
g.ax_heatmap.hlines(y=len(table2.index), xmin=-0.5, xmax=len(table2.columns), linewidth=5, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table2.index), linewidth=5, color='black')
g.ax_heatmap.vlines(x=len(table2.columns), ymin=0, ymax=len(table2.index), linewidth=5, color='black')
g.cax.set_position([3.95,1.6, 0.07, 0.8]) 
g.cax.set_title("RMSD",loc="center",pad=20,fontsize=35)
g.cax.tick_params(labelsize=35)
g.ax_row_colors.tick_params(labelsize=30)
g.ax_col_colors.tick_params(labelsize=30)

#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(2.2,-0.2),fontsize=35,title_fontsize=35, bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=5,title='Gprotein family/GPCR class')


g.savefig(path+fil+'htmp_classA.svg',bbox_inches='tight',dpi=300)




#violin plots full
dg=df2[(df2['gfam1']=='Gs') & (df2['gfam2']=='Gs')]
di=df2[(df2['gfam1']=='Gio') & (df2['gfam2']=='Gio')]
dq=df2[(df2['gfam1']=='Gq11') & (df2['gfam2']=='Gq11')]

dgi=pd.concat([dg,di])
dgi.to_csv("./plots/violin_Gs_Gi_all"+fil+".tsv",index=None,header=True,sep='\t')

sns.set(style="ticks",rc={'axes.facecolor':'white'},font_scale=2)

fig , ax = plt.subplots(figsize=(10,10))
ax = sns.violinplot(x=dgi["gfam1"],y=dgi["rmsd"],palette=lut, dodge=False)
ax.yaxis.set_major_locator(plt.MaxNLocator(10))


pairs=[('Gs', 'Gio')]
stat,pval=stats.ranksums(dg["rmsd"], di["rmsd"])
#stat1,pval1=stats.ranksums(dg["rmsd"], dq["rmsd"])
#stat2,pval2=stats.ranksums(di["rmsd"], dq["rmsd"])
test_short_name = 'Ranksum'
pvalues = [pval]
#ax.set_ylim(y_lims[0], 1.25 * y_lims[1])

ax.set_ylim([-3, 20])  
ax.set_xlabel("G-protein", fontsize=50)
ax.set_ylabel("RMSD",fontsize=50)
ax.tick_params(labelsize=50)
#plt.yticks(np.arange(-1, dgi['rmsd'].max(),10))
annot = Annotator(ax, pairs,plot='violinplot', data=dgi, x='gfam1', y='rmsd')
(annot.configure(test=None, test_short_name=test_short_name, text_format='full', loc='outside').set_pvalues(pvalues=pvalues).annotate())
#y_lims = ax.get_ylim()
plt.savefig(path+fil+'vlp_all.svg',bbox_inches='tight')

with open(path+fil+'_stats.txt', 'w') as f:
    print('Gprotein',fs,'Median',fs,'std',fs,'sem',file=f)
    print ('Gs',fs,np.median(dg["rmsd"]),fs, np.std(dg["rmsd"]),fs, stats.sem(dg["rmsd"]),file=f)
    print ('Gio',fs,np.median(di["rmsd"]), np.std(di["rmsd"]), stats.sem(di["rmsd"]),file=f)
    #print ('Gq11',fs,np.median(dq["rmsd"]), np.std(dq["rmsd"]), stats.sem(dq["rmsd"]),file=f)



#class A
#violin plots 
dg=df3[(df3['gfam1']=='Gs') & (df3['gfam2']=='Gs')]
di=df3[(df3['gfam1']=='Gio') & (df3['gfam2']=='Gio')]
dq=df3[(df3['gfam1']=='Gq11') & (df3['gfam2']=='Gq11')]

dgi=pd.concat([dg,di])
dgi.to_csv("./plots/violin_Gs_GiclassA_"+fil+".tsv",index=None,header=True,sep='\t')

sns.set(style="ticks",rc={'axes.facecolor':'white'},font_scale=2)

fig , ax = plt.subplots(figsize=(10,10))
ax = sns.violinplot(x=dgi["gfam1"],y=dgi["rmsd"],palette=lut, dodge=False)
ax.yaxis.set_major_locator(plt.MaxNLocator(10))


pairs=[('Gs', 'Gio')]
stat,pval=stats.ranksums(dg["rmsd"], di["rmsd"])
#stat1,pval1=stats.ranksums(dg["rmsd"], dq["rmsd"])
#stat2,pval2=stats.ranksums(di["rmsd"], dq["rmsd"])
test_short_name = 'Ranksum'
pvalues = [pval]
#ax.set_ylim(y_lims[0], 1.25 * y_lims[1])

ax.set_ylim([-3,20])
ax.set_xlabel("G-protein", fontsize=50)
ax.set_ylabel("RMSD",fontsize=50)
ax.tick_params(labelsize=50)
#plt.yticks(np.arange(-1, dgi['rmsd'].max(),10))
annot = Annotator(ax, pairs,plot='violinplot', data=dgi, x='gfam1', y='rmsd')
(annot.configure(test=None, test_short_name=test_short_name, text_format='full', loc='outside').set_pvalues(pvalues=pvalues).annotate())
#y_lims = ax.get_ylim()
plt.savefig(path+fil+'vlp_classA.svg',bbox_inches='tight')

with open(path+fil+'_stats_classA.txt', 'w') as f:
    print('Gprotein',fs,'Median',fs,'std',fs,'sem',file=f)
    print ('Gs',fs,np.median(dg["rmsd"]),fs, np.std(dg["rmsd"]),fs, stats.sem(dg["rmsd"]),file=f)
    print ('Gio',fs,np.median(di["rmsd"]), np.std(di["rmsd"]), stats.sem(di["rmsd"]),file=f)
    #print ('Gq11',fs,np.median(dq["rmsd"]), np.std(dq["rmsd"]), stats.sem(dq["rmsd"]),file=f)

