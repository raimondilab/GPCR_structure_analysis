
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


############### Calculate and plot log odds ratio of Gs vs Gio structure contacts for paired GPCR-Gprotein contacts




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

#pair logodds
alld=df1
gprotcount=df1.groupby(['gfam'])['GPCR'].nunique().reset_index()
alld['unique_gpcr_contact'] = alld.groupby(['comb'])['GPCR'].transform('nunique')
s=alld.groupby(['comb','unique_gpcr_contact']).size().reset_index(name='Total_count')
alld1=alld.merge(gprotcount,left_on=['gfam'],right_on=['gfam'],how="left")
alld1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
alld1['gpcr_coupled_contact'] = alld1.groupby(['comb','gfam'])['GPCR_x'].transform('nunique')
alld1['gpcr_coupled_nocontact']=alld1['coupled_gpcr']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_contact']=alld1['unique_gpcr_contact']-alld1['gpcr_coupled_contact']
alld1['gpcr_uncoupled_nocontact']=alld1['GPCR_x'].nunique()-alld1['coupled_gpcr']-alld1['gpcr_uncoupled_contact']
alld1.to_csv('./heatmaps/logs_pair/log_odds_pairs_all_result.tsv',index=True,header=True,sep='\t')


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
fin2=fin1.loc[abs(fin1['ln(odds_ratio)']) > 1]
fin3=fin1.loc[abs(fin1['ln(odds_ratio)']) > 2]

t1 = pd.pivot_table(fin1, values='ln(odds_ratio)', index=['Gfam'],columns=['comb'])
t1=t1.T
t2 = pd.pivot_table(fin2, values='ln(odds_ratio)', index=['Gfam'],columns=['comb'])
t2=t2.T
t3 = pd.pivot_table(fin3, values='ln(odds_ratio)', index=['Gfam'],columns=['comb'])
t3=t3.T
#lst=t1.columns.tolist()
#t1[lst] = scaler.fit_transform(t1[lst])
tf=t1
tf1=t2
tf2=t3
tf.to_csv('./heatmaps/logs_pair/log_odds_pairs_all.tsv',index=True,header=True,sep='\t')
tf1.to_csv('./heatmaps/logs_pair/log_odds_pairs_abs1.tsv',index=True,header=True,sep='\t')
tf2.to_csv('./heatmaps/logs_pair/log_odds_pairs_abs2.tsv',index=True,header=True,sep='\t')



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
plt.ylabel('Pair positions(%)', fontsize=40)
plt.yticks(np.arange(0,110,10))
plt.legend(bbox_to_anchor=(0.147,1),ncol=2,fontsize=35,title='log_odds_ratio',title_fontsize=40)
plt.savefig('./heatmaps/logs_pair/barplot_count_gprot_all.svg',dpi=600,bbox_inches='tight',format='svg')


####01
gpcrs=df1.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrs = gpcrs[(gpcrs['Total_count']/len(set(df1['GPCR'])) > 0.1)]
df1gc=gpcrs.merge(df1,left_on=['comb'],right_on=['comb'],how="left")
lr1g=df1gc["comb"].drop_duplicates().tolist()
tf01 = tf.reindex(lr1g)
tf011 = tf1.reindex(lr1g).dropna(how='all')
tf012 = tf2.reindex(lr1g).dropna(how='all')

tf01.to_csv('./heatmaps/logs_pair/log_odds_pairs_01.tsv',index=True,header=True,sep='\t')
tf011.to_csv('./heatmaps/logs_pair/log_odds_pairs_01_abs1.tsv',index=True,header=True,sep='\t')
tf012.to_csv('./heatmaps/logs_pair/log_odds_pairs_01_abs2.tsv',index=True,header=True,sep='\t')

#classA
classA=df1.loc[df1['Class'] == 'classA']
gpcrsA=classA.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrsAA = gpcrsA[(gpcrsA['Total_count']/len(set(classA['GPCR'])) > 0.1)]
df1gcA=gpcrsAA.merge(classA,left_on=['comb'],right_on=['comb'],how="left")
l1A=classA["comb"].drop_duplicates().tolist()
lr1A=df1gcA["comb"].drop_duplicates().tolist()
tfA = tf.reindex(l1A)
tfA01 = tf.reindex(lr1A)
tfA.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_all_classA.tsv',index=True,header=True,sep='\t')
tfA01.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_01_classA.tsv',index=True,header=True,sep='\t')
tfA1 = tf1.reindex(l1A).dropna(how='all')
tfA2 = tf2.reindex(l1A).dropna(how='all')
tfA1.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_all_classA_abs1.tsv',index=True,header=True,sep='\t')
tfA2.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_all_classA_abs2.tsv',index=True,header=True,sep='\t')
tfA011 = tf1.reindex(lr1A).dropna(how='all')
tfA012 = tf2.reindex(lr1A).dropna(how='all')
tfA011.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_01_classA_abs1.tsv',index=True,header=True,sep='\t')
tfA012.to_csv('./heatmaps/logs_pair/log_odds_pair_pos_01_classA_abs2.tsv',index=True,header=True,sep='\t')






############################################# LOG PLOT
###gpcr
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
fig.savefig('./heatmaps/logs_pair/log_odds_pair.svg',dpi=600,bbox_inches='tight',format='svg')
