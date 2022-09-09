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
df=pd.read_csv("all_gpcrs.tsv", comment="#", sep="\t")
clas=pd.read_csv("classification.txt", comment="#", sep="\t")
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
#gs=df[df.gfam =='Gs']
#gq11=df[df.gfam =='Gq11']
gs['val']=1
gio['val']=1
#gq11['val']=3
#df1=pd.concat([gs,gio,gq11],ignore_index=True)
df1=pd.concat([gs,gio],ignore_index=True)

#gprot logodds
alld=df1
gprotcount=df1.groupby(['gfam'])['GPCR'].nunique().reset_index()
alld['GPCR']=np.where(((alld['GPCR']=='HTR4') & (alld['gfam']=='Gio')),'HTR4-Gio',alld['GPCR'])
alld['GPCR']=np.where(((alld['GPCR']=='CCKAR') & (alld['gfam']=='Gio')),'CCKAR-Gio',alld['GPCR'])
alld['GPCR']=np.where(((alld['GPCR']=='GCGR') & (alld['gfam']=='Gio')),'GCGR-Gio',alld['GPCR'])
alld['GPCR']=np.where(((alld['GPCR']=='CCKBR') & (alld['gfam']=='Gio')),'CCKBR-Gio',alld['GPCR'])
alld['GPCR']=np.where(((alld['GPCR']=='GPR139') & (alld['gfam']=='Gio')),'GPR139-Gio',alld['GPCR'])
alld['unique_gpcr_contact'] = alld.groupby(['comb'])['GPCR'].transform('nunique')
s=alld.groupby(['comb','unique_gpcr_contact']).size().reset_index(name='Total_count')
alld1=alld.merge(gprotcount,left_on=['gfam'],right_on=['gfam'],how="left")
alld1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
alld1['gpcr_coupled_contact'] = alld1.groupby(['comb','gfam'])['GPCR_x'].transform('nunique')
alld1['gpcr_coupled_nocontact']=alld1['coupled_gpcr']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_contact']=alld1['unique_gpcr_contact']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_nocontact']=alld1['GPCR_x'].nunique()-alld1['coupled_gpcr']-alld1['gpcr_uncoupled_contact']
alld1.to_csv('log_odds_gprot_pair_result.tsv',index=True,header=True,sep='\t')

dic={}
for index, row in alld1.iterrows():
    dic[index]=[row['comb'],row['gfam'],np.asarray([[row['gpcr_coupled_contact'], row['gpcr_coupled_nocontact']],[row['gpcr_uncoupled_contact'],row['gpcr_uncoupled_nocontact']]])]
dic1={} 
for key,value in dic.items():
    dic1[key]=[value[0],value[1],sm.stats.Table2x2(value[2],shift_zeros=True)]
dic2={}
for key,value in dic1.items():
    dic2[key] =[value[0],value[1],value[2].oddsratio,value[2].log_oddsratio,sm.stats.Table2x2.log_oddsratio_pvalue(value[2])]  
fin1=pd.DataFrame.from_dict(dic2,orient='index',columns=['comb','Gfam','odds_ratio', 'ln(odds_ratio)', 'ln(odds_ratio)_pval'])
fin1=fin1.drop_duplicates()
fin2=fin1.loc[abs(fin1['ln(odds_ratio)']) > 2]
t1 = pd.pivot_table(fin1, values='ln(odds_ratio)', index=['Gfam'],columns=['comb'])
t1=t1.T
t2 = pd.pivot_table(fin2, values='ln(odds_ratio)', index=['Gfam'],columns=['comb'])
t2=t2.T
#lst=t1.columns.tolist()
#t1[lst] = scaler.fit_transform(t1[lst])
tf=t1
tf1=t2
tf.to_csv('log_odds_pair_pos_new.tsv',index=True,header=True,sep='\t')
tf1.to_csv('log_odds_pair_pos_abs2.tsv',index=True,header=True,sep='\t')


#bar count of position pairs per coupled gpcr that meet conditions above results gfamily strucutural element
br1=t1.T.unstack().reset_index().dropna()

br1['activity']=pd.cut(br1[0],[-7,-1,0,1,7], labels=["worse", "bad","good","better"])
br1=br1[['Gfam','activity']]
br1['<-1'] = br1.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('worse',case=True)].count())
br1['-1-0'] = br1.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('bad', case=True)].count())
br1['0-1'] = br1.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('good', case=True)].count())
br1['>1'] = br1.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('better', case=True)].count())
br1['Gfam']=br1['Gfam'].astype(str).str.strip()
br1['sum']=br1['<-1'] +br1['-1-0']+br1['0-1']+br1['>1']
br1[['<-1','-1-0','0-1','>1']] = (br1[['<-1','-1-0','0-1','>1']].div(br1['sum'], axis=0)*100)
count=br1[['Gfam','<-1','-1-0','0-1','>1']].drop_duplicates()

count.set_index(['Gfam'],inplace=True)
sns.set_style('ticks')
fig = plt.figure(figsize=(10,10))
fig.set_facecolor('white')
count.reindex(["Gs", "Gio"]).plot.bar(stacked=True, figsize=(10, 10), cmap='RdYlGn', edgecolor='None',fontsize=40,width=0.85)
plt.xlabel('Gprotein', fontsize=40)
plt.ylabel('Gprotein positions(%)', fontsize=40)
plt.yticks(np.arange(0,110,10))
plt.legend(bbox_to_anchor=(0.147,1),ncol=2,fontsize=35,title='log_odds_ratio',title_fontsize=40)
plt.savefig('barplot_count_pair_all_unscaled.svg',dpi=600,bbox_inches='tight',format='svg')


##########
##########
df1= df1.merge(clas, left_on=['GPCR'], right_on=['GPCR'],how="left")
df1['GPCR']=np.where(((df1['GPCR']=='HTR4') & (df1['gfam']=='Gio')),'HTR4-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKAR') & (df1['gfam']=='Gio')),'CCKAR-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='GCGR') & (df1['gfam']=='Gio')),'GCGR-Gio',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKAR') & (df1['gfam']=='Gs')),'CCKAR-Gs',df1['GPCR'])
df1['GPCR']=np.where(((df1['GPCR']=='CCKBR') & (df1['gfam']=='Gio')),'CCKBR-Gio',df1['GPCR'])
gpcrs=df1.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gprots=df1.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.1)]
df1gc=gpcrs.merge(df1,left_on=['comb'],right_on=['comb'],how="left")

sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1.Structural = df1.Structural.astype("category")
df1.Structural.cat.set_categories(sorter, inplace=True)
df1.sort_values(["Structural","comb"],inplace=True)
lr1=df1["comb"].drop_duplicates().tolist()


#comb
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr1)

t3=t2.reset_index()
comb=t3[['comb']]
comb01 = comb.merge(df1gc,left_on=['comb'],right_on=['comb'],how="left")
comb01=comb01.dropna()
df1gc=comb01
#02
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1gc.Structural = df1gc.Structural.astype("category")
df1gc.Structural.cat.set_categories(sorter, inplace=True)
df1gc.sort_values(["Structural","comb"],inplace=True)
lr1g=df1gc["comb"].drop_duplicates().tolist()

#comb01
table12g = pd.pivot_table(df1gc ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
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
clas=pd.read_csv("classification.txt", comment="#", sep="\t")
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
fam=pd.read_csv("targets_and_families.csv", comment="#", sep=",")
fam= gpcr.merge(fam, left_on=['GPCRs'], right_on=['GPCR'],how="left")
fam=fam.drop(fam.columns[[1,2,3,4]], axis=1).drop_duplicates()
fam=fam.set_index('GPCR_x')
fam = fam.pop("Family")
fam=fam.astype('object')
lut2 = dict(zip(set(fam), sns.husl_palette(len(set(fam)),h=.5)))
rs2= fam.map(lut2)
coup=pd.read_csv("coup_assay.txt", sep="\t",index_col=0)

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
#combplot
cl=ListedColormap(['olive'])
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=5.5)
g=sns.clustermap(table12g,mask=table12g==mask,row_colors=rs3, tree_kws=dict(linewidths=2, colors=(0, 0, 0)),method='ward',figsize=(33,33),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.1), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.cax.set_visible(False)
g.gs.update(top=2.6,right=2)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
ax2 = divider.append_axes("top", size="8000%", pad=0)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g.ax_heatmap.hlines(y=len(table12g.index), xmin=-0.5, xmax=len(table12g.columns), linewidth=20, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
g.ax_heatmap.vlines(x=len(table12g.columns), ymin=0, ymax=len(table12g.index), linewidth=20, color='black')
t2 = t2.reindex(lr1g)
t2.to_csv('log_odds_pair_pos_abs2_01.tsv',index=True,header=True,sep='\t')

t2.plot.bar(stacked=True, color=['b','r'], edgecolor='black',fontsize=15,width=0.85,legend=False,ax=ax2)
ax2.set_xticklabels([])
ax2.set_xlabel('')
ax2.set_facecolor('white')
ax2.tick_params(axis=u'both', which=u'both',length=1, direction="inout")
ax2.hlines(y=0, xmin=-0.5, xmax=len(table12g.columns)-0.5, linewidth=4, color='black')
ax2.hlines(y=-5, xmin=-0.5, xmax=len(table12g.columns)-0.5, linewidth=4, color='black')
ax2.hlines(y=5, xmin=-0.5, xmax=len(table12g.columns)-0.5, linewidth=4, color='black')
ax2.vlines(x=-0.5, ymin=-5, ymax=5, linewidth=4, color='black')
ax2.vlines(x=len(table12g.columns)-0.52, ymin=-5, ymax=5, linewidth=4, color='black')
ax2.set_xlim(-0.5,len(table12g.columns)-0.5)
plt.ylabel('log_odds_ratio',fontsize=65,labelpad=30)
plt.tick_params(axis='y', which='both', labelsize=55)
plt.yticks(np.arange(-5,5,1))
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g.ax_row_dendrogram.legend(loc="best",ncol=3, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(2.4,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
    leg1= g.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.9,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text().split('-')[0] in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text().split('-')[0]])


g.savefig('Gs_vs_Gio_receptor_gprot_cluster_pair_01_lodd2.svg',dpi=600,bbox_inches='tight',format='svg')


####classA
df1=df1.loc[df1['class'] == 'classA']
df1=df1[(df.BW != 'ICL1') & (df1.BW != 'Cterm') & (df.BW != 'ICL2') & (df1.BW != 'ICL3') & (df1.BW != 'H8')]

gpcrs=df1.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.1)]
df1gc=gpcrs.merge(df1,left_on=['comb'],right_on=['comb'],how="left")


sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1.Structural = df1.Structural.astype("category")
df1.Structural.cat.set_categories(sorter, inplace=True)
df1.sort_values(["Structural","comb"],inplace=True)
lr1=df1["comb"].drop_duplicates().tolist()

#gprot
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr1)
mask =0
table1=table1.dropna(axis=1, how='all')
t3=tf1.reset_index()
comb=t3[['comb']]
comb01 = comb.merge(df1gc,left_on=['comb'],right_on=['comb'],how="left")
comb01=comb01.dropna()
df1gc=comb01
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
df1gc.Structural = df1gc.Structural.astype("category")
df1gc.Structural.cat.set_categories(sorter, inplace=True)
df1gc.sort_values(["Structural","comb"],inplace=True)
lr1g=df1gc["comb"].drop_duplicates().tolist()

table1g = pd.pivot_table(df1gc ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
table1g = table1g.reindex(columns=lr1g)
table1g=table1g.dropna(axis=1, how='all')

mask =0
# create the customized color map gfamily
ds=df1[['GPCR', 'gfam']].drop_duplicates()
ds.set_index('GPCR',inplace=True)
ds = ds.pop("gfam")
ds=ds.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}

rs = ds.map(lut)
# create the customized color map gpcr class
clas=pd.read_csv("classification.txt", comment="#", sep="\t")
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
fam=pd.read_csv("targets_and_families.csv", comment="#", sep=",")
fam= gpcr.merge(fam, left_on=['GPCRs'], right_on=['GPCR'],how="left")
fam=fam.drop(fam.columns[[1,2,3,4]], axis=1).drop_duplicates()
fam=fam.set_index('GPCR_x')
fam = fam.pop("Family")
fam=fam.astype('object')
lut2 = dict(zip(set(fam), sns.husl_palette(len(set(fam)),h=.5)))
rs2= fam.map(lut2)
coup=pd.read_csv("coup_assay.txt", sep="\t",index_col=0)

rs3 = pd.concat([rs,rs1,rs2],axis=1)
rs3 = pd.concat([rs3,coup],axis=1)
lut.update(lut1)
#coup={'G1213':'#800080','Only GtoPdb':'#cfccc6','Only TGF or GEMTA':'#FF8C00','Not coupled':'#808080'}
coup={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}
lut.update(coup)





#combplot
cl=ListedColormap(['olive'])
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
sns.set(font_scale=5.5)
g=sns.clustermap(table1g,mask=table1g==mask,row_colors=rs3, tree_kws=dict(linewidths=2, colors=(0, 0, 0)),method='ward',figsize=(33,33),dendrogram_ratio=(.05, .1),colors_ratio=(0.01,0.1), col_cluster=False,cmap=cl,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver')
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =65)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 65) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.cax.set_visible(False)
g.gs.update(top=2.6,right=2)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
ax2 = divider.append_axes("top", size="8000%", pad=0)
g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g.ax_heatmap.hlines(y=len(table1g.index), xmin=-0.5, xmax=len(table1g.columns), linewidth=20, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
g.ax_heatmap.vlines(x=len(table1g.columns), ymin=0, ymax=len(table1g.index), linewidth=20, color='black')
t2 = tf1.reindex(lr1g)
t2.to_csv('log_odds_pair_pos_abs2_01_classA.tsv',index=True,header=True,sep='\t')

t2.plot.bar(stacked=True, color=['b','r'], edgecolor='black',fontsize=15,width=0.85,legend=False,ax=ax2)
ax2.set_xticklabels([])
ax2.set_xlabel('')
ax2.set_facecolor('white')
ax2.tick_params(axis=u'both', which=u'both',length=1, direction="inout")
ax2.hlines(y=0, xmin=-0.5, xmax=len(table1g.columns)-0.5, linewidth=4, color='black')
ax2.hlines(y=-5, xmin=-0.5, xmax=len(table1g.columns)-0.5, linewidth=4, color='black')
ax2.hlines(y=5, xmin=-0.5, xmax=len(table1g.columns)-0.5, linewidth=4, color='black')
ax2.vlines(x=-0.5, ymin=-5, ymax=5, linewidth=4, color='black')
ax2.vlines(x=len(table1g.columns)-0.52, ymin=-5, ymax=5, linewidth=4, color='black')
ax2.set_xlim(-0.5,len(table1g.columns)-0.5)
plt.ylabel('log_odds_ratio',fontsize=65,labelpad=30)
plt.tick_params(axis='y', which='both', labelsize=55)
plt.yticks(np.arange(-5,5,1))
#gfamily and class
handles = [Patch(facecolor=lut[name]) for name in lut]
leg2=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,0.07), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=2,fontsize=70,title='Gprotein family(UCM)/GPCR class',title_fontsize=70,columnspacing=0.1,handletextpad =0.1)
#gpcr family
for label in fam.unique():
    g.ax_row_dendrogram.bar(0, 0, color=lut2[label], label=label, linewidth=0)
leg1=g.ax_row_dendrogram.legend(loc="best",ncol=3, title='GPCR family',bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=60,bbox_to_anchor=(2.4,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct element
for label,v in ann_colors.items():
    g.ax_col_dendrogram.bar(0, 0, color=v, label=label, linewidth=0)
    leg1= g.ax_col_dendrogram.legend(title='Structural element',loc="best",ncol=2, bbox_transform=plt.gcf().transFigure,frameon=False,fontsize=70,bbox_to_anchor=(0.9,0.07),columnspacing=0.2,title_fontsize=70,handletextpad =0.2)
#struct leg
for xtic in g.ax_heatmap.get_xmajorticklabels():
    if xtic.get_text().split('-')[0] in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text().split('-')[0]])

g.savefig('Gs_vs_Gio_receptor_gprot_cluster_pair_01_lodd2_classA.svg',dpi=600,bbox_inches='tight',format='svg')
