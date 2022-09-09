#!/data/SW/anaconda3/envs/myenv/bin/python

from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.preprocessing import MaxAbsScaler
from mpl_toolkits.axes_grid1 import make_axes_locatable


#Gs vs Gi pairs
df=pd.read_csv("./use_file/all_gpcrs_cons.tsv", comment="#", sep="\t")
clas=pd.read_csv("./use_file/classification.txt", comment="#", sep="\t")
df.rename(columns = {'GPCR_x':'GPCR'}, inplace = True)
df[['G', 'Gprot_struct','num']] = df['CGN'].str.split('.', expand=True)
df['Gprot_pos']=df['Gprot_struct']+"."+df["num"]
df=df.drop(df.columns[[4,7]], axis=1)
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein"].map(equiv)
df['comb']=df["BW"]+"-"+df["Gprot_pos"]
df['num']=df['num'].astype(float)
ann_colors={"Nterm":"crimson",
            "Cterm":"deepskyblue",
                "TM1":"orange",
               "TM2":"yellowgreen",
               "TM3":"olive",
               "TM4":"green",
               "TM5":"lime",
               "TM6":"gold",
               "TM7":"darkkhaki",
               "H8":"black",
               "ECL1":"darkorchid",
               "ECL2":"gray",
               "ECL3":"teal",
               "ICL1":"navy",
               "ICL2":"tan",
               "ICL3":"fuchsia"}


d1 = pd.DataFrame([[key,value] for key,value in ann_colors.items()],
     columns=["Structural","hex"])

df = df.merge(d1,left_on=['Structural'],right_on=['Structural'],how="left")
t=df.iloc[:,[4,11]] 
t.columns = ['Structural', 'hex']
mydict = dict(zip(t.Structural, t.hex))
#by the receptor
gs=df[df.gfam =='Gs']
gio=df[df.gfam =='Gio']
gq11=df[df.gfam =='Gq11']
g1213=df[df.gfam =='G1213']
gs['val']=1
gio['val']=1
gq11['val']=1
g1213['val']=1

df1=pd.concat([gs,gio,gq11,g1213],ignore_index=True)
#df1=pd.concat([gs,gio],ignore_index=True)

##########
df1= df1.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
df1['GPCR']=np.where(((df1['GPCR']=='HTR4') & (df1['gfam']=='Gio')),'HTR4-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKAR') & (df1['gfam']=='Gio')),'CCKAR-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='GCGR') & (df1['gfam']=='Gio')),'GCGR-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKAR') & (df1['gfam']=='Gs')),'CCKAR-Gs',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKBR') & (df1['gfam']=='Gio')),'CCKBR-Gio',df1['GPCR'])
gpcrs=df1.groupby(['BW'])['GPCR'].nunique().reset_index(name='Total_count')
gprots=df1.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.2)]
gprots = gprots[(gprots['Total_count']/len(set(df1['GPCR'])) > 0.2)]
df1gc=gpcrs.merge(df1,left_on=['BW'],right_on=['BW'],how="left")
df1gp=gprots.merge(df1,left_on=['Gprot_pos'],right_on=['Gprot_pos'],how="left")

sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1.Structural = df1.Structural.astype("category")
df1.Structural.cat.set_categories(sorter, inplace=True)
df1.sort_values(["Structural","BW"],inplace=True)
lr1=df1["BW"].drop_duplicates().tolist()
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
df1.Gprot_struct = df1.Gprot_struct.astype("category")
df1.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1=df1[['Gprot_struct','Gprot_pos','num']].drop_duplicates()
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2=ds1["Gprot_pos"].drop_duplicates().tolist()
sorter2=['Gs','Gio','Gq11','G1213']
#gprot
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr2)
#gpcr
table12 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12= table12.reindex(columns=lr1)
mask =0

#02
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1gc.Structural = df1gc.Structural.astype("category")
df1gc.Structural.cat.set_categories(sorter, inplace=True)
df1gc.sort_values(["Structural","BW"],inplace=True)
lr1g=df1gc["BW"].drop_duplicates().tolist()
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
df1gp.Gprot_struct = df1gp.Gprot_struct.astype("category")
df1gp.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1=df1gp[['Gprot_struct','Gprot_pos','num']].drop_duplicates()
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2g=ds1["Gprot_pos"].drop_duplicates().tolist()
sorter2=['Gs','Gio','Gq11','G1213']
#gprot
table1g = pd.pivot_table(df1gp ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
table1g = table1g.reindex(columns=lr2g)
#gpcr
table12g = pd.pivot_table(df1gc ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12g= table12g.reindex(columns=lr1g)
mask =0

# create the customized color map gfamily
ds=df1[['GPCR', 'gfam']].drop_duplicates()
ds.set_index('GPCR',inplace=True)
ds = ds.pop("gfam")
ds=ds.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}
rs = ds.map(lut)
# create the customized color map gpcr class
clas=pd.read_csv("./use_file/classification.txt", comment="#", sep="\t")
gpcr=ds.reset_index()
gpcr[['GPCRs', 'lin']] = gpcr['GPCR'].str.split('-', 1, expand=True)
clas= gpcr.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
clas=clas.drop(clas.columns[[1,2,3]], axis=1).drop_duplicates()
clas=clas.set_index('GPCR')
clas = clas.pop("class")
clas=clas.astype('object')
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}
rs1 = clas.map(lut1)
# create the customized color map gpcr family
fam=pd.read_csv("./use_file/targets_and_families.csv", comment="#", sep=",")
fam= gpcr.merge(fam, left_on=['GPCRs'], right_on=['GPCR'],how="left")
fam=fam.drop(fam.columns[[1,2,3,4]], axis=1).drop_duplicates()
fam=fam.set_index('GPCR_x')
fam = fam.pop("Family")
fam=fam.astype('object')
lut2 = dict(zip(set(fam), sns.husl_palette(len(set(fam)),h=.5)))
rs2= fam.map(lut2)
coup=pd.read_csv("./use_file/coup_assay.txt", sep="\t",index_col=0)

rs3 = pd.concat([rs,rs1,rs2],axis=1)
rs3 = pd.concat([rs3,coup],axis=1)

lut.update(lut1)
#coup={'G1213':'#800080','Only GtoPdb':'#cfccc6','Only TGF or GEMTA':'#FF8C00','Not coupled':'#808080'}
coup={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}

lut.update(coup) 

#gprot ann
gprot_ann={'Helix':'deepskyblue','Sheet':'crimson','Loop':'darkorchid'}
gprotp=df1['Gprot_pos'].to_frame().drop_duplicates()
gprotp1 = gprotp[~gprotp['Gprot_pos'].str.isupper()]
gprotp1['col']='darkorchid'
gprotp2 = gprotp[gprotp['Gprot_pos'].str.isupper() & gprotp['Gprot_pos'].str.startswith('H')]
gprotp3 = gprotp[gprotp['Gprot_pos'].str.isupper() & gprotp['Gprot_pos'].str.startswith('S')]
gprotp2['col']='deepskyblue'
gprotp3['col']='crimson'
gp=pd.concat([gprotp1,gprotp2,gprotp3])
gp_ann = dict(zip(gp.Gprot_pos, gp.col))

#gprot plot
#cl=ListedColormap(['r','b','#008000','#800080'])
cl=ListedColormap(['olive'])

sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=8)
g=sns.clustermap(table1,mask=table1==mask,row_colors=rs3, tree_kws=dict(linewidths=2, colors=(0, 0, 0)),method='ward',figsize=(63,30),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.1), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.cax.set_visible(False)
g.gs.update(top=2.6,right=2)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1.columns), linewidth=20, color='black')
g.ax_heatmap.hlines(y=len(table1.index), xmin=-0.5, xmax=len(table1.columns), linewidth=20, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1.index), linewidth=20, color='black')
g.ax_heatmap.vlines(x=len(table1.columns), ymin=0, ymax=len(table1.index), linewidth=20, color='black')
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g.ax_row_dendrogram.legend(loc="best",ncol=4, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(1.91,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in gprot_ann.items():
    g.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=1, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.8,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
        xtic.set_color(gp_ann[xtic.get_text()])

g.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gprot_cluster_new.svg',dpi=600,bbox_inches='tight',format='svg')

#### gpcr
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=8)
g1=sns.clustermap(table12,mask=table12==mask,row_colors=rs3,method='ward',tree_kws=dict(linewidths=2, colors=(0, 0, 0)),figsize=(66,30),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.01), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g1.ax_heatmap.set_xticklabels(g1.ax_heatmap.get_xmajorticklabels(), fontsize =63)
g1.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g1.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g1.ax_heatmap.set_yticklabels(g1.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g1.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g1.cax.set_visible(False)
g1.gs.update(top=2.5,right=2)
g1.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table12.columns), linewidth=20, color='black')
g1.ax_heatmap.hlines(y=len(table12.index), xmin=-0.5, xmax=len(table12.columns), linewidth=20, color='black')
g1.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table12.index), linewidth=20, color='black')
g1.ax_heatmap.vlines(x=len(table12.columns), ymin=0, ymax=len(table12.index), linewidth=20, color='black')
a= g1.ax_heatmap
a.set_xlabel("GPCR Residue number")
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g1.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g1.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g1.ax_row_dendrogram.legend(loc="best",ncol=4, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(1.9,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g1.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
    leg1= g1.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.8,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g1.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])

g1.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gpcr_cluster_new.svg',dpi=600,bbox_inches='tight',format='svg')


#gprot plot
#cl=ListedColormap(['r','b','#008000','#800080'])
cl=ListedColormap(['olive'])

sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=6)
g3=sns.clustermap(table1g,mask=table1g==mask,method='ward',row_colors=rs3,figsize=(20,35),tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.05, .1),colors_ratio=(0.025,0.05), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g3.ax_heatmap.set_xlabel('Position', fontsize=70)
g3.ax_heatmap.set_ylabel('GPCR', fontsize=70)
g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g3.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g3.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g3.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g3.cax.set_visible(False)
g3.gs.update(top=2,right=2)
g3.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g3.ax_heatmap.hlines(y=len(table1g.index), xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g3.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
g3.ax_heatmap.vlines(x=len(table1g.columns), ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g3.ax_heatmap.legend(handles, lut,bbox_to_anchor=(1.2,0.02), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=3,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g3.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
    leg1=g3.ax_row_dendrogram.legend(loc="best",ncol=1, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(3.45,1.8),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
for label,v in gprot_ann.items():
    g3.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
    leg1= g3.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=1, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.7,0.02),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g3.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
        xtic.set_color(gp_ann[xtic.get_text()])

g3.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gprot_cluster_02.svg',dpi=600,bbox_inches='tight',format='svg')

#### gpcr
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=6)
g4=sns.clustermap(table12g,mask=table12g==mask,method='ward',row_colors=rs3,figsize=(20,35),tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.05, .1),colors_ratio=(0.025,0.05), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g4.ax_heatmap.set_xlabel('Position', fontsize=70)
g4.ax_heatmap.set_ylabel('GPCR', fontsize=70)
g4.ax_heatmap.set_xticklabels(g4.ax_heatmap.get_xmajorticklabels(), fontsize =62)
g4.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g4.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g4.ax_heatmap.set_yticklabels(g4.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g4.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g4.cax.set_visible(False)
g4.gs.update(top=2,right=2)
g4.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g4.ax_heatmap.hlines(y=len(table12g.index), xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g4.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
g4.ax_heatmap.vlines(x=len(table12g.columns), ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
a= g4.ax_heatmap
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g4.ax_heatmap.legend(handles, lut,bbox_to_anchor=(1.2,0.02), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=3,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g4.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
    leg1=g4.ax_row_dendrogram.legend(loc="best",ncol=1, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(3.3,1.8),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g4.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g4.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.75,0.02),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g4.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
g4.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gpcr_cluster_02.svg',dpi=600,bbox_inches='tight',format='svg')


####classA
#gprot
#df1= df1.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
df1=df1.loc[df1['class'] == 'classA']
gpcrs=df1.groupby(['BW'])['GPCR'].nunique().reset_index(name='Total_count')
gprots=df1.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.2)]
gprots = gprots[(gprots['Total_count']/len(set(df1['GPCR'])) > 0.2)]
df1gc=gpcrs.merge(df1,left_on=['BW'],right_on=['BW'],how="left")
df1gp=gprots.merge(df1,left_on=['Gprot_pos'],right_on=['Gprot_pos'],how="left")


sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1.Structural = df1.Structural.astype("category")
df1.Structural.cat.set_categories(sorter, inplace=True)
df1.sort_values(["Structural","BW"],inplace=True)
lr1=df1["BW"].drop_duplicates().tolist()
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
df1.Gprot_struct = df1.Gprot_struct.astype("category")
df1.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1=df1[['Gprot_struct','Gprot_pos','num']].drop_duplicates()
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2=ds1["Gprot_pos"].drop_duplicates().tolist()
sorter2=['Gs','Gio','Gq11','G1213']
#gprot
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr2)
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr2)
mask =0
table1=table1.dropna(axis=1, how='all')
#gpcr
table12 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12= table12.reindex(columns=lr1)
table12=table12.dropna(axis=1, how='all')
mask =0
table12=table12.drop(['ICL1','Cterm','H8','ICL2','ICL3'], axis=1)
lr1.remove('ICL1')
lr1.remove('ICL2')
lr1.remove('ICL3')
lr1.remove('H8')
lr1.remove('Cterm')
#gprot

###02

sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1gc.Structural = df1gc.Structural.astype("category")
df1gc.Structural.cat.set_categories(sorter, inplace=True)
df1gc.sort_values(["Structural","BW"],inplace=True)
lr1g=df1gc["BW"].drop_duplicates().tolist()
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
df1gp.Gprot_struct = df1gp.Gprot_struct.astype("category")
df1gp.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1=df1gp[['Gprot_struct','Gprot_pos','num']].drop_duplicates()
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2g=ds1["Gprot_pos"].drop_duplicates().tolist()
sorter2=['Gs','Gio','Gq11','G1213']


table1g = pd.pivot_table(df1gp ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
table1g = table1g.reindex(columns=lr2g)
table1g=table1g.dropna(axis=1, how='all')
#gpcr
#df1gc= df1gc.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
#02

table12g = pd.pivot_table(df1gc ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12g= table12g.reindex(columns=lr1g)
table12g=table12g.dropna(axis=1, how='all')
table12g=table12g.drop(['ICL3'], axis=1)
lr1g.remove('ICL3')

mask =0

# create the customized color map gfamily
ds=df1[['GPCR', 'gfam']].drop_duplicates()
ds.set_index('GPCR',inplace=True)
ds = ds.pop("gfam")
ds=ds.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}

rs = ds.map(lut)
# create the customized color map gpcr class
clas=pd.read_csv("./use_file/classification.txt", comment="#", sep="\t")
gpcr=ds.reset_index()
gpcr[['GPCRs', 'lin']] = gpcr['GPCR'].str.split('-', 1, expand=True)
clas= gpcr.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
clas=clas.drop(clas.columns[[1,2,3]], axis=1).drop_duplicates()
clas=clas.set_index('GPCR')
clas = clas.pop("class")
clas=clas.astype('object')
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}
lut1 = {'classA':'#44daeb'}

rs1 = clas.map(lut1)
# create the customized color map gpcr family
fam=pd.read_csv("./use_file/targets_and_families.csv", comment="#", sep=",")
fam= gpcr.merge(fam, left_on=['GPCRs'], right_on=['GPCR'],how="left")
fam=fam.drop(fam.columns[[1,2,3,4]], axis=1).drop_duplicates()
fam=fam.set_index('GPCR_x')
fam = fam.pop("Family")
fam=fam.astype('object')
lut2 = dict(zip(set(fam), sns.husl_palette(len(set(fam)),h=.5)))
rs2= fam.map(lut2)
coup=pd.read_csv("./use_file/coup_assay.txt", sep="\t",index_col=0)

rs3 = pd.concat([rs,rs1,rs2],axis=1)
rs3 = pd.concat([rs3,coup],axis=1)

lut.update(lut1)
#coup={'G1213':'#800080','Only GtoPdb':'#cfccc6','Only TGF or GEMTA':'#FF8C00','Not coupled':'#808080'}
coup={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}

lut.update(coup) 


#gprot plot
#cl=ListedColormap(['r','b','#008000','#800080'])
cl=ListedColormap(['olive'])

sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=9)
g=sns.clustermap(table1,mask=table1==mask,row_colors=rs3, tree_kws=dict(linewidths=2, colors=(0, 0, 0)),method='ward',figsize=(63,30),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.1), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =70)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 68) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.cax.set_visible(False)
g.gs.update(top=2.6,right=2)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1.columns), linewidth=20, color='black')
g.ax_heatmap.hlines(y=len(table1.index), xmin=-0.5, xmax=len(table1.columns), linewidth=20, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1.index), linewidth=20, color='black')
g.ax_heatmap.vlines(x=len(table1.columns), ymin=0, ymax=len(table1.index), linewidth=20, color='black')
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g.ax_row_dendrogram.legend(loc="best",ncol=4, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.99,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in gprot_ann.items():
    g.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=1, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.8,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
        xtic.set_color(gp_ann[xtic.get_text()])

g.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gprot_cluster_classA.svg',dpi=600,bbox_inches='tight',format='svg')

#### gpcr
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=9)
g1=sns.clustermap(table12,mask=table12==mask,row_colors=rs3,method='ward',tree_kws=dict(linewidths=2, colors=(0, 0, 0)),figsize=(66,30),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.01), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g1.ax_heatmap.set_xticklabels(g1.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g1.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g1.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g1.ax_heatmap.set_yticklabels(g1.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g1.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g1.cax.set_visible(False)
g1.gs.update(top=2.5,right=2)
g1.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table12.columns), linewidth=20, color='black')
g1.ax_heatmap.hlines(y=len(table12.index), xmin=-0.5, xmax=len(table12.columns), linewidth=20, color='black')
g1.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table12.index), linewidth=20, color='black')
g1.ax_heatmap.vlines(x=len(table12.columns), ymin=0, ymax=len(table12.index), linewidth=20, color='black')
a= g1.ax_heatmap
a.set_xlabel("GPCR Residue number")
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g1.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g1.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g1.ax_row_dendrogram.legend(loc="best",ncol=4, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.99,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g1.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g1.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.8,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g1.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])

g1.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gpcr_cluster_classA.svg',dpi=600,bbox_inches='tight',format='svg')


#gprot plot
#cl=ListedColormap(['r','b','#008000','#800080'])
cl=ListedColormap(['olive'])

sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=6)
g3=sns.clustermap(table1g,mask=table1g==mask,method='ward',row_colors=rs3,figsize=(20,35),tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.05, .1),colors_ratio=(0.025,0.05), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g3.ax_heatmap.set_xlabel('Position', fontsize=70)
g3.ax_heatmap.set_ylabel('GPCR', fontsize=70)
g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g3.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g3.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g3.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g3.cax.set_visible(False)
g3.gs.update(top=2,right=2)
g3.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g3.ax_heatmap.hlines(y=len(table1g.index), xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g3.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
g3.ax_heatmap.vlines(x=len(table1g.columns), ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g3.ax_heatmap.legend(handles, lut,bbox_to_anchor=(1.2,0.02), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=3,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g3.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
    leg1=g3.ax_row_dendrogram.legend(loc="best",ncol=1, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(3.45,1.8),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
for label,v in gprot_ann.items():
    g3.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g3.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=1, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.7,0.02),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g3.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
        xtic.set_color(gp_ann[xtic.get_text()])

g3.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gprot_cluster_classA_02.svg',dpi=600,bbox_inches='tight',format='svg')

#### gpcr
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=6)
g2=sns.clustermap(table12g,mask=table12g==mask,method='ward',row_colors=rs3,figsize=(20,35),tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.05, .1),colors_ratio=(0.025,0.05), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g2.ax_heatmap.set_xlabel('Position', fontsize=70)
g2.ax_heatmap.set_ylabel('GPCR', fontsize=70)
g2.ax_heatmap.set_xticklabels(g2.ax_heatmap.get_xmajorticklabels(), fontsize =62)
g2.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g2.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g2.ax_heatmap.set_yticklabels(g2.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g2.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g2.cax.set_visible(False)
g2.gs.update(top=2,right=2)
g2.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g2.ax_heatmap.hlines(y=len(table12g.index), xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g2.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
g2.ax_heatmap.vlines(x=len(table12g.columns), ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
a= g2.ax_heatmap
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g2.ax_heatmap.legend(handles, lut,bbox_to_anchor=(1.2,0.02), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=3,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g2.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g2.ax_row_dendrogram.legend(loc="best",ncol=1, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(3.3,1.8),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g2.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
leg1= g2.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(1.75,0.02),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g2.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
g2.savefig('./heatmaps/plots1/Gs_vs_Gio_receptor_gpcr_cluster_classA_02.svg',dpi=600,bbox_inches='tight',format='svg')

