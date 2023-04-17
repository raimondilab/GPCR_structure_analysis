
from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
from pylab import *
import sys




############### Calculate and plot log odds ratio of Gs vs Gio structure contacts separatley for GPCR and Gprotein portion

#Gs vs Gi pairs
df=pd.read_csv("./use_file/all_gpcrs.tsv", comment="#", sep="\t")
infile = open('./use_file/GPCR_class.pickle', 'rb')
clas = pickle.load(infile)
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein"].map(equiv)
df['comb']=df["BW"]+"-"+df["Gprot_pos"]
df["Class"] = df["GPCR"].map(clas)

gs=df[df.gfam =='Gs']
gio=df[df.gfam =='Gio']
#gs=df[df.gfam =='Gs']
#gq11=df[df.gfam =='Gq11']
gs['val']=1
gio['val']=1
#gq11['val']=3
#df1=pd.concat([gs,gio,gq11],ignore_index=True)
df1=pd.concat([gs,gio],ignore_index=True)
###gpcrs with multiple coupling

dupl = df1[['GPCR', 'gfam']].drop_duplicates()
dupl_mask = dupl['GPCR'].duplicated()
gpcr_dupl = dupl[dupl_mask]
for index,row in gpcr_dupl.iterrows():
    df1['GPCR']=np.where(((df1['GPCR']==row['GPCR']) & (df1['gfam']==row['gfam'])),row['GPCR']+"-"+row['gfam'],df1['GPCR'])

#gprot logodds
alld=df1
gprotcount=df1.groupby(['gfam'])['GPCR'].nunique().reset_index()
alld['unique_gpcr_contact'] = alld.groupby(['Gprot_pos'])['GPCR'].transform('nunique')
s=alld.groupby(['Gprot_pos','unique_gpcr_contact']).size().reset_index(name='Total_count')
alld1=alld.merge(gprotcount,left_on=['gfam'],right_on=['gfam'],how="left")
alld1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
alld1['gpcr_coupled_contact'] = alld1.groupby(['Gprot_pos','gfam'])['GPCR_x'].transform('nunique')
alld1['gpcr_coupled_nocontact']=alld1['coupled_gpcr']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_contact']=alld1['unique_gpcr_contact']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_nocontact']=alld1['GPCR_x'].nunique()-alld1['coupled_gpcr']-alld1['gpcr_uncoupled_contact']
alld1.to_csv('./heatmaps/logs/log_odds_gprot_all_result.tsv',index=True,header=True,sep='\t')

dic={}
for index, row in alld1.iterrows():
    dic[index]=[row['Gprot_pos'],row['gfam'],np.asarray([[row['gpcr_coupled_contact'], row['gpcr_coupled_nocontact']],[row['gpcr_uncoupled_contact'],row['gpcr_uncoupled_nocontact']]])]
dic1={} 
for key,value in dic.items():
    dic1[key]=[value[0],value[1],sm.stats.Table2x2(value[2],shift_zeros=True)]
dic2={}
for key,value in dic1.items():
    dic2[key] =[value[0],value[1],value[2].oddsratio,value[2].log_oddsratio,sm.stats.Table2x2.log_oddsratio_pvalue(value[2])]  
fin1=pd.DataFrame.from_dict(dic2,orient='index',columns=['Gprot_pos','Gfam','odds_ratio', 'ln(odds_ratio)', 'ln(odds_ratio)_pval'])
fin1=fin1.drop_duplicates()
t1 = pd.pivot_table(fin1, values='ln(odds_ratio)', index=['Gfam'],columns=['Gprot_pos'])
t1=t1.T
tf=t1
tf1=t1
tf.to_csv('./heatmaps/logs/log_odds_gprot_pos_all.tsv',index=True,header=True,sep='\t')

#gpcr logodds
ald=df1
gprtcount=df1.groupby(['gfam'])['GPCR'].nunique().reset_index()
ald['unique_gpcr_contact'] = alld.groupby(['Gprot_pos'])['GPCR'].transform('nunique')
ald['unique_gpcr_contact'] = ald.groupby(['BW'])['GPCR'].transform('nunique')
s=ald.groupby(['Gprot_pos','unique_gpcr_contact']).size().reset_index(name='Total_count')
ald1=ald.merge(gprtcount,left_on=['gfam'],right_on=['gfam'],how="left")
ald1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
ald1['gpcr_coupled_contact'] = ald1.groupby(['BW','gfam'])['GPCR_x'].transform('nunique')
ald1['gpcr_coupled_nocontact']=ald1['coupled_gpcr']-ald1['gpcr_coupled_contact']
ald1['gpcr_uncoupled_contact']=ald1['unique_gpcr_contact']-ald1['gpcr_coupled_contact']
ald1['gpcr_uncoupled_nocontact']=ald1['GPCR_x'].nunique()-ald1['coupled_gpcr']-ald1['gpcr_uncoupled_contact']
ald1.to_csv('./heatmaps/logs/log_odds_gpcr_all_result.tsv',index=True,header=True,sep='\t')
di={}
for index, row in ald1.iterrows():
    di[index]=[row['BW'],row['gfam'],np.asarray([[row['gpcr_coupled_contact'], row['gpcr_coupled_nocontact']],[row['gpcr_uncoupled_contact'],row['gpcr_uncoupled_nocontact']]])]
di1={} 
for key,value in di.items():
    di1[key]=[value[0],value[1],sm.stats.Table2x2(value[2],shift_zeros=True)]
di2={}
for key,value in di1.items():
    di2[key] =[value[0],value[1],value[2].oddsratio,value[2].log_oddsratio,sm.stats.Table2x2.log_oddsratio_pvalue(value[2])]  
fi1=pd.DataFrame.from_dict(di2,orient='index',columns=['BW','Gfam','odds_ratio', 'ln(odds_ratio)', 'ln(odds_ratio)_pval'])
fi1=fi1.drop_duplicates()
tg = pd.pivot_table(fi1, values='ln(odds_ratio)', index=['Gfam'],columns=['BW'])
tg=tg.T
lst1=tg.columns.tolist()
tg1=tg
tg.to_csv('./heatmaps/logs/log_odds_gpcr_pos_all.tsv',index=True,header=True,sep='\t')

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
plt.savefig('./heatmaps/logs/barplot_count_gprot_all.svg',dpi=600,bbox_inches='tight',format='svg')
#gpcr
br2=tg.T.unstack().reset_index().dropna()
br2['activity']=pd.cut(br2[0],[-7,-1,0,1,7], labels=["worse", "bad","good","better"])

br2=br2[['Gfam','activity']]
br2['<-1'] = br2.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('worse',case=True)].count())
br2['-1-0'] = br2.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('bad', case=True)].count())
br2['0-1'] = br2.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('good', case=True)].count())
br2['>1'] = br2.groupby(['Gfam'])['activity'].transform(lambda x: x[x.str.contains('better', case=True)].count())
br2['Gfam']=br2['Gfam'].astype(str).str.strip()
br2['sum']=br2['<-1'] +br2['-1-0']+br2['0-1']+br2['>1']
br2[['<-1','-1-0','0-1','>1']] = (br2[['<-1','-1-0','0-1','>1']].div(br2['sum'], axis=0)*100)
count1=br2[['Gfam','<-1','-1-0','0-1','>1']].drop_duplicates()
count1.set_index(['Gfam'],inplace=True)

sns.set_style('ticks')
fig = plt.figure(figsize=(10,10))
fig.set_facecolor('white')

count1.reindex(["Gs", "Gio"]).plot.bar(stacked=True, figsize=(10, 10), cmap='RdYlGn', edgecolor='None',fontsize=40,width=0.85)
plt.xlabel('Gprotein', fontsize=40)
plt.ylabel('GPCR positions(%)', fontsize=40)
plt.yticks(np.arange(0,110,10))
plt.legend(bbox_to_anchor=(0.15, 1),fontsize=35,ncol=2,title='log_odds_ratio',title_fontsize=40)
plt.savefig('./heatmaps/logs/barplot_count1_gpcr_all_unscaled.svg',dpi=600,bbox_inches='tight',format='svg')


####02
gpcrs=df1.groupby(['BW'])['GPCR'].nunique().reset_index(name='Total_count')
gprots=df1.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.2)]
gprots = gprots[(gprots['Total_count']/len(set(df1['GPCR'])) > 0.2)]
df1gc=gpcrs.merge(df1,left_on=['BW'],right_on=['BW'],how="left")
df1gp=gprots.merge(df1,left_on=['Gprot_pos'],right_on=['Gprot_pos'],how="left")
lr1g=df1gc["BW"].drop_duplicates().tolist()
lr2g=df1gp["Gprot_pos"].drop_duplicates().tolist()
tf02 = tf.reindex(lr2g)
tg02 = tg.reindex(lr1g)
tg02.to_csv('./heatmaps/logs/log_odds_gpcr_pos_02.tsv',index=True,header=True,sep='\t')
tf02.to_csv('./heatmaps/logs/log_odds_gprot_pos_02.tsv',index=True,header=True,sep='\t')

#classA
classA=df1.loc[df1['Class'] == 'classA']
gpcrsA=classA.groupby(['BW'])['GPCR'].nunique().reset_index(name='Total_count')
gprotsA=classA.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')

gpcrsAA = gpcrsA[(gpcrsA['Total_count']/len(set(classA['GPCR'])) > 0.2)]
gprotsAA = gprotsA[(gprotsA['Total_count']/len(set(classA['GPCR'])) > 0.2)]
df1gcA=gpcrsAA.merge(classA,left_on=['BW'],right_on=['BW'],how="left")
df1gpA=gprotsAA.merge(classA,left_on=['Gprot_pos'],right_on=['Gprot_pos'],how="left")
l1A=classA["BW"].drop_duplicates().tolist()
l2A=classA["Gprot_pos"].drop_duplicates().tolist()
lr1A=df1gcA["BW"].drop_duplicates().tolist()
lr2A=df1gpA["Gprot_pos"].drop_duplicates().tolist()
tfA = tf.reindex(l2A)
tgA = tg.reindex(l1A)
tfA02 = tf.reindex(lr2A)
tgA02 = tg.reindex(lr1A)
tgA.to_csv('./heatmaps/logs/log_odds_gpcr_pos_all_classA.tsv',index=True,header=True,sep='\t')
tfA.to_csv('./heatmaps/logs/log_odds_gprot_pos_all_classA.tsv',index=True,header=True,sep='\t')
tgA02.to_csv('./heatmaps/logs/log_odds_gpcr_pos_02_classA.tsv',index=True,header=True,sep='\t')
tfA02.to_csv('./heatmaps/logs/log_odds_gprot_pos_02_classA.tsv',index=True,header=True,sep='\t')
############################################# LOG PLOT
###gpcr
sns.set(style="ticks",rc={'axes.facecolor':'white','figure.figsize': (8, 5)})
fig, (ax, ax1) = plt.subplots(nrows=2, sharex=True)
plt.subplots_adjust(wspace=0, hspace=0)
subplot(2,1,1)
ax=sns.kdeplot(data=tg, x='Gio')
y= ax.lines[-1].get_ydata()
x=ax.lines[-1].get_xdata()

#ax.set_xlim(-6.5, 6.5)
med=tg['Gio'].median()
ax.vlines(med, 0,y.max(), linestyle='--', alpha=1,color='black')
plt.text(0.46, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(left=False)
ax.axes.yaxis.set_visible(False)
ax.set_xticks(np.arange(-7,7,1))
plt.xticks(rotation=90)
plt.fill_between(x,0,y, where = x<=-1, color='#a50026')
plt.fill_between(x,0,y, where = x>=1, color='#006837')
plt.fill_between(x,0,y, where = (x>=0) & (x<=1), color='#b7e075')
plt.fill_between(x,y, where = (x>=-1) & (x<=0), color='#fdbf6f')
plt.text(0.1, 0.5, 'Gio', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=25)
plt.gca().invert_xaxis()

subplot(2,1,2)
ax1=sns.kdeplot(data=tg, x='Gs')
ax1.axes.yaxis.set_ticklabels([])
ax1.axes.yaxis.set_visible(False)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(left=False)
ax1.set_xlim(-7, 7)
ax1.set_xticks(np.arange(-7,7,0.5))
plt.xticks(rotation=90)
med=tg['Gs'].median()
y1= ax1.lines[-1].get_ydata()
x1=ax1.lines[-1].get_xdata()
ax1.vlines(med, 0,y1.max(), linestyle='--', alpha=1,color='black')
plt.fill_between(x1,0,y1, where = x1<-1, color='#a50026')
plt.fill_between(x1,0,y1, where = x1>1, color='#006837')
plt.fill_between(x1,0,y1, where = (x1>=0) & (x1<1), color='#b7e075')
plt.fill_between(x1,0,y1, where = (x1>-1) & (x1<=0), color='#fdbf6f')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.xlabel("logodds",fontsize=25,rotation=180)
plt.text(0.1, 0.5, 'Gs', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=25)
plt.text(0.45, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=15)
fig.savefig('./heatmaps/logs/log_odds_gpcr.svg',dpi=600,bbox_inches='tight',format='svg')
stat,pval=stats.ranksums(tg["Gs"], tg["Gio"])
print(pval)
#gprot
sns.set(style="ticks",rc={'axes.facecolor':'white','figure.figsize': (8, 5)})
fig, (ax, ax1) = plt.subplots(nrows=2, sharex=True)
plt.subplots_adjust(wspace=0, hspace=0)
subplot(2,1,1)
ax=sns.kdeplot(data=tf, x='Gio')
y= ax.lines[-1].get_ydata()
x=ax.lines[-1].get_xdata()

#ax.set_xlim(-6.5, 6.5)
med=tf['Gio'].median()
ax.vlines(med, 0,y.max(), linestyle='--', alpha=1,color='black')
plt.text(0.46, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=15)
ax.axes.yaxis.set_ticklabels([])
ax.tick_params(left=False)
ax.axes.yaxis.set_visible(False)
ax.set_xticks(np.arange(-7,7,1))
plt.xticks(rotation=90)
plt.fill_between(x,0,y, where = x<=-1, color='#a50026')
plt.fill_between(x,0,y, where = x>=1, color='#006837')
plt.fill_between(x,0,y, where = (x>=0) & (x<=1), color='#b7e075')
plt.fill_between(x,y, where = (x>=-1) & (x<=0), color='#fdbf6f')
plt.text(0.1, 0.5, 'Gio', horizontalalignment='center',
     verticalalignment='center', transform=ax.transAxes,rotation=90,fontsize=25)
plt.gca().invert_xaxis()

subplot(2,1,2)
ax1=sns.kdeplot(data=tf, x='Gs')
ax1.axes.yaxis.set_ticklabels([])
ax1.axes.yaxis.set_visible(False)
ax1.tick_params(axis='x', labelsize=20)
ax1.tick_params(left=False)
ax1.set_xlim(-7, 7)
ax1.set_xticks(np.arange(-7,7,0.5))
plt.xticks(rotation=90)
med=tf['Gs'].median()
y1= ax1.lines[-1].get_ydata()
x1=ax1.lines[-1].get_xdata()
ax1.vlines(med, 0,y1.max(), linestyle='--', alpha=1,color='black')
plt.fill_between(x1,0,y1, where = x1<-1, color='#a50026')
plt.fill_between(x1,0,y1, where = x1>1, color='#006837')
plt.fill_between(x1,0,y1, where = (x1>=0) & (x1<1), color='#b7e075')
plt.fill_between(x1,0,y1, where = (x1>-1) & (x1<=0), color='#fdbf6f')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.xlabel("logodds",fontsize=25,rotation=180)
plt.text(0.1, 0.5, 'Gs', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=25)
plt.text(0.4, 0.5, 'median', horizontalalignment='center',
     verticalalignment='center', transform=ax1.transAxes,rotation=90,fontsize=15)
fig.savefig('./heatmaps/logs/'+'log_odds_gprot.svg',dpi=600,bbox_inches='tight',format='svg')
stat1,pval1=stats.ranksums(tf["Gs"], tf["Gio"])
print(pval1)
################################################ mapping log odds to gprotein and gpcr positions
'''

gprot=pd.read_csv("human_sequence_table_FR_new.txt", comment="#", sep="\t",index_col=0)
gprot1=gprot.unstack().reset_index()
gprot1.columns=['Gprotein_id','CGN','Pos2']
gprot1['Pos2'] = gprot1['Pos2'].str.replace(r'\D', '')
gprot1['Pos2'] .replace('', np.nan, inplace=True)
gprot1=gprot1.dropna()
gpcr=pd.read_csv("output_gpcrdb.tsv", dtype=object, sep="\t")
df=pd.read_csv("GPCR_structs_clean.tsv", sep="\t", dtype=object)
d=pd.read_csv("uniprot_pdb.tsv", sep="\t", dtype=object)
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein_Gene_name"].map(equiv)
gpcr_a=d.merge(gpcr,left_on=['UNIPROT_AC','UNIPROT_RES_NUM'],right_on=['Uniprot','Pos1'],how="left")
gpcr_a=gpcr_a.dropna()
gfam=df[['PDB_ID','gfam']]
#gpcrs
loggpcr=pd.read_csv("log_odds_gpcr_pos_new.tsv", dtype=object, sep="\t")
loggpcr=loggpcr.set_index('BW')
loggpcr=loggpcr.unstack().reset_index().dropna()
loggpcr.columns=['gfam','BW','logodds']
gpcr_a1=gpcr_a.merge(gfam,left_on=['#PDBID'],right_on=['PDB_ID'],how="left")
gpcrs=gpcr_a1[['PDB_ID','CHAIN','PDB_RES_NUM','UNIPROT_RES_NUM','BW','Structural','Uniprot','GPCR','gfam']]
gpcrs=gpcrs.merge(loggpcr,left_on=['gfam','BW'],right_on=['gfam','BW'],how="left")
gpcrs=gpcrs.dropna()
gpcrs.to_csv('gpcr_loggodds_mapping_cons.txt',index=False,header=True,sep='\t')

#gprot
gfam1=df[['PDB_ID','Gprotein_Uniprot_ID','Gprotein_Uniprot_AC','gfam']]
gprot2=gprot1.merge(gfam1,left_on=['Gprotein_id'],right_on=['Gprotein_Uniprot_ID'],how="left")
gprot2=gprot2.dropna()
gprot_a=d.merge(gprot2,left_on=['#PDBID','UNIPROT_AC','UNIPROT_RES_NUM'],right_on=['PDB_ID','Gprotein_Uniprot_AC','Pos2'],how="left")
gprot_a=gprot_a.dropna()
loggprot=pd.read_csv("log_odds_gprot_pos_new.tsv", dtype=object, sep="\t")
loggprot=loggprot.set_index('Gprot_pos')
loggprot=loggprot.unstack().reset_index().dropna()
loggprot.columns=['gfam','CGN','logodds']
gprot_a1=gprot_a[['PDB_ID','CHAIN','PDB_RES_NUM','UNIPROT_RES_NUM','CGN','UNIPROT_AC','Gprotein_id','gfam']]
gprot_a1['CGN']=gprot_a1["CGN"].str.replace("G.","")
gprot_a1=gprot_a1.merge(loggprot,left_on=['gfam','CGN'],right_on=['gfam','CGN'],how="left")
gprot_a1=gprot_a1.dropna()
gpcrs.to_csv('gprot_loggodds_mapping_cons.txt',index=False,header=True,sep='\t')

'''
