#!/data/SW/anaconda3/envs/myenv/bin/python


import statsmodels.api as sm
from sklearn.preprocessing import MaxAbsScaler
from matplotlib.patches import Patch

import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import pandas as pd
import numpy as np
fil=sys.argv[1]
df=pd.read_csv("./network/files/gio_contact_"+fil+".tsv", comment="#", sep="\t")
df = df[(df['unique_gpcr_contact']/len(set(df['GPCR'])) > 0.25)]
df1=pd.read_csv("./network/files/gs_contact_"+fil+".tsv", comment="#", sep="\t")
df1 = df1[(df1['unique_gpcr_contact']/len(set(df1['GPCR'])) > 0.25)]
#df2=pd.read_csv("G1213_contact.tsv", comment="#", sep="\t")
df3=pd.read_csv("./network/files/gq11_contact_"+fil+".tsv", comment="#", sep="\t")
df3 = df3[(df3['unique_gpcr_contact']/len(set(df3['GPCR'])) > 0.25)]

d_all=pd.concat([df,df1,df3]).drop_duplicates(keep=False)
d_all = d_all[d_all["BW"].str.contains("Cterm") ==False]

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

d_all = d_all.merge(d1,left_on=['Structural'],right_on=['Structural'],how="left")
t=d_all.iloc[:,[4,17]] 
t.columns = ['Structural', 'hex']
mydict = dict(zip(t.Structural, t.hex))


gprot_ann={'Helix':'deepskyblue','Sheet':'crimson','Loop':'darkorchid'}
gprotp=d_all['Gprot_pos'].to_frame().drop_duplicates()
gprotp1 = gprotp[~gprotp['Gprot_pos'].str.isupper()]
gprotp1['col']='darkorchid'
gprotp2 = gprotp[gprotp['Gprot_pos'].str.isupper() & gprotp['Gprot_pos'].str.startswith('H')]
gprotp3 = gprotp[gprotp['Gprot_pos'].str.isupper() & gprotp['Gprot_pos'].str.startswith('S')]
gprotp2['col']='deepskyblue'
gprotp3['col']='crimson'
gp=pd.concat([gprotp1,gprotp2,gprotp3])
gp_ann = dict(zip(gp.Gprot_pos, gp.col))


gss=d_all[d_all.gfam =='Gs']
gio=d_all[d_all.gfam =='Gio']
#g1213=d_all[d_all.gfam =='G1213']
gq11=d_all[d_all.gfam =='Gq11']
#gio
pos=d_all[['BW','Gprot_pos','Structural','Gprot_struct']].drop_duplicates()
gio = pos.merge(gio,on=['BW','Gprot_pos','Structural','Gprot_struct'],how="left")
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
gio.Gprot_struct = gio.Gprot_struct.astype("category")
gio.Gprot_struct.cat.set_categories(sorter1, inplace=True)
gio.sort_values(["Gprot_struct","Gprot_pos"],inplace=True)
ls2=gio["Gprot_pos"].drop_duplicates().tolist()
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
gio.Structural = gio.Structural.astype("category")
gio.Structural.cat.set_categories(sorter, inplace=True)
gio.sort_values(["Structural","BW"],inplace=True)
ls1=gio["BW"].drop_duplicates().tolist()
gio['Gprot_pos'] = pd.Categorical(gio['Gprot_pos'],
                                   categories=ls2,
                                   ordered=True)


table = pd.pivot_table(gio, values='fraction_pair', index=['Gprot_pos'],columns=['BW'])
table= table.reindex(ls2)
table = table.reindex(columns=ls1)
#gs
gss = pos.merge(gss,on=['BW','Gprot_pos','Structural','Gprot_struct'],how="left")
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
gss.Gprot_struct = gss.Gprot_struct.astype("category")
gss.Gprot_struct.cat.set_categories(sorter1, inplace=True)
gss.sort_values(["Gprot_struct","Gprot_pos"],inplace=True)
ls3=gss["Gprot_pos"].drop_duplicates().tolist()
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
gss.Structural = gss.Structural.astype("category")
gss.Structural.cat.set_categories(sorter, inplace=True)
gss.sort_values(["Structural","BW"],inplace=True)
ls4=gss["BW"].drop_duplicates().tolist()
gss['Gprot_pos'] = pd.Categorical(gss['Gprot_pos'],
                                   categories=ls3,
                                   ordered=True)

table1 = pd.pivot_table(gss, values='fraction_pair', index=['Gprot_pos'],columns=['BW'])
table1= table1.reindex(ls3)
table1 = table1.reindex(columns=ls4)

'''
#G1213
g1213 = pos.merge(g1213,on=['BW','CGN','Structural','Struct'],how="left")
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
g1213.Struct = g1213.Struct.astype("category")
g1213.Struct.cat.set_categories(sorter1, inplace=True)
g1213.sort_values(["Struct","CGN"],inplace=True)
ls5=g1213["CGN"].drop_duplicates().tolist()
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
g1213.Structural = g1213.Structural.astype("category")
g1213.Structural.cat.set_categories(sorter, inplace=True)
g1213.sort_values(["Structural","BW"],inplace=True)
ls6=g1213["BW"].drop_duplicates().tolist()
g1213['CGN'] = pd.Categorical(g1213['CGN'],
                                   categories=ls5,
                                   ordered=True)

table2 = pd.pivot_table(g1213, values='fraction_pair', index=['CGN'],columns=['BW'])
table2= table2.reindex(ls5)
table2 = table2.reindex(columns=ls6)
'''
#gq11
gq11 = pos.merge(gq11,on=['BW','Gprot_pos','Structural','Gprot_struct'],how="left")
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
gq11.Gprot_struct = gq11.Gprot_struct.astype("category")
gq11.Gprot_struct.cat.set_categories(sorter1, inplace=True)
gq11.sort_values(["Gprot_struct","Gprot_pos"],inplace=True)
ls5=gq11["Gprot_pos"].drop_duplicates().tolist()
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
gq11.Structural = gq11.Structural.astype("category")
gq11.Structural.cat.set_categories(sorter, inplace=True)
gq11.sort_values(["Structural","BW"],inplace=True)
ls6=gq11["BW"].drop_duplicates().tolist()
gq11['Gprot_pos'] = pd.Categorical(gq11['Gprot_pos'],
                                   categories=ls5,
                                   ordered=True)
table3= pd.pivot_table(gq11, values='fraction_pair', index=['Gprot_pos'],columns=['BW'])
table3= table3.reindex(ls5)
table3 = table3.reindex(columns=ls6)
#####
alld=d_all.iloc[:, [0,1,2,3,4,5,6,7,8]]

gpcrcount=d_all.groupby(['gfam'])['GPCR'].nunique().reset_index()
alld['unique_gpcr_contact'] = alld.groupby(['Gprot_pos'])['GPCR'].transform('nunique')
s=alld.groupby(['Gprot_pos','unique_gpcr_contact']).size().reset_index(name='Total_count')
alld1=alld.merge(gpcrcount,left_on=['gfam'],right_on=['gfam'],how="left")
alld1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
alld1['gpcr_coupled_contact'] = alld1.groupby(['Gprot_pos','gfam'])['GPCR_x'].transform('nunique')
alld1['gpcr_coupled_nocontact']=alld1['coupled_gpcr']-alld1['gpcr_coupled_contact']
alld1['gpcr_unoupled_contact']=alld1['unique_gpcr_contact']-alld1['gpcr_coupled_contact']
alld1['gpcr_unoupled_nocontact']=alld1['GPCR_x'].nunique()-alld1['coupled_gpcr']-alld1['gpcr_unoupled_contact']
dic={}
for index, row in alld1.iterrows():
    dic[index]=[row['Gprot_pos'],row['gfam'],np.asarray([[row['gpcr_coupled_contact'], row['gpcr_coupled_nocontact']],[row['gpcr_unoupled_contact'],row['gpcr_unoupled_nocontact']]])]
dic1={} 
for key,value in dic.items():
    dic1[key]=[value[0],value[1],sm.stats.Table2x2(value[2],shift_zeros=True)]
dic2={}
for key,value in dic1.items():
    dic2[key] =[value[0],value[1],value[2].oddsratio,value[2].log_oddsratio,sm.stats.Table2x2.log_oddsratio_pvalue(value[2])]  
fin1=pd.DataFrame.from_dict(dic2,orient='index',columns=['Gprot_pos','Gfam','odds_ratio', 'ln(odds_ratio)', 'ln(odds_ratio)_pval'])
fin1=fin1.drop_duplicates()
scaler = MaxAbsScaler()
t1 = pd.pivot_table(fin1, values='ln(odds_ratio)', index=['Gfam'],columns=['Gprot_pos'])
t1=t1.T
lst=t1.columns.tolist()
t1[lst] = scaler.fit_transform(t1[lst])
tf=t1
tf.to_csv('./heatmaps/files/log_odds_gprot_pos_025_'+fil+'.tsv',index=True,header=True,sep='\t')

alld=d_all.iloc[:, [0,1,2,3,4,5,6,7,8]]
alld['unique_gpcr_contact'] = alld.groupby(['BW'])['GPCR'].transform('nunique')
s=alld.groupby(['BW','unique_gpcr_contact']).size().reset_index(name='Total_count')
alld1=alld.merge(gpcrcount,left_on=['gfam'],right_on=['gfam'],how="left")
alld1.rename(columns = {'GPCR_y':'coupled_gpcr'}, inplace = True)
alld1['gpcr_coupled_contact'] = alld1.groupby(['BW','gfam'])['GPCR_x'].transform('nunique')
alld1['gpcr_coupled_nocontact']=alld1['coupled_gpcr']-alld1['gpcr_coupled_contact']
alld1['gpcr_unoupled_contact']=alld1['unique_gpcr_contact']-alld1['gpcr_coupled_contact']
alld1['gpcr_unoupled_nocontact']=alld1['GPCR_x'].nunique()-alld1['coupled_gpcr']-alld1['gpcr_unoupled_contact']
dic={}
for index, row in alld1.iterrows():
    dic[index]=[row['BW'],row['gfam'],np.asarray([[row['gpcr_coupled_contact'], row['gpcr_coupled_nocontact']],[row['gpcr_unoupled_contact'],row['gpcr_unoupled_nocontact']]])]
dic1={} 
for key,value in dic.items():
    dic1[key]=[value[0],value[1],sm.stats.Table2x2(value[2],shift_zeros=True)]
dic2={}
for key,value in dic1.items():
    dic2[key] =[value[0],value[1],value[2].oddsratio,value[2].log_oddsratio,sm.stats.Table2x2.log_oddsratio_pvalue(value[2])]  
fin1=pd.DataFrame.from_dict(dic2,orient='index',columns=['BW','Gfam','odds_ratio', 'ln(odds_ratio)', 'ln(odds_ratio)_pval'])
fin1=fin1.drop_duplicates()
scaler = MaxAbsScaler()
t1 = pd.pivot_table(fin1, values='ln(odds_ratio)', index=['Gfam'],columns=['BW'])
t1=t1.T
lst=t1.columns.tolist()
t1[lst] = scaler.fit_transform(t1[lst])
tg=t1
tg.to_csv('./heatmaps/files/log_odds_gprot_pos_025_'+fil+'.tsv',index=True,header=True,sep='\t')
tg = tg.reindex(ls1)
tg=tg.reset_index().fillna(0)
tf = tf.reindex(ls2)
tf=tf.reset_index().fillna(0)
#Gio
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'},font_scale=5)
g=sns.clustermap(table,cmap='Blues',row_cluster=False,col_cluster=False,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver',vmin=0,vmax=1,figsize=(20,20),cbar_kws=dict(use_gridspec=False,pad=10,ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]))
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =50)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 50) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.gs.update(top=1.5,right=1.5)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
ax= divider.append_axes("top", size="12000%", pad=0)
ax2= divider2.append_axes("right", size="12000%", pad=0)
tg.plot.bar(x='BW',y='Gio', color=['b'],ax=ax ,edgecolor='black',width=0.9,legend=False)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_facecolor('white')
ax.set_ylabel('log_odds_ratio',fontsize=50,labelpad=30)
ax.tick_params(axis='y', which='both', labelsize=42)
ax.set_yticks(np.arange(-1,1.2,0.2))
ax.set_ylim(-1.2,1)
tf.plot.barh(x='Gprot_pos',y='Gio',ax=ax2, color=['b'] ,edgecolor='black',fontsize=15,width=0.9,legend=False)
ax2.set_yticklabels([])
ax2.set_ylabel('')
ax2.set_facecolor('white')
ax2.set_xlabel('log_odds_ratio',fontsize=50,labelpad=20)
ax2.tick_params(axis='x', which='both', labelsize=42,rotation=90)
ax2.set_xlim(-1.2,1)
ax2.set_xticks(np.arange(-1,1.2,0.2))
ax2.invert_xaxis()
g.ax_heatmap.hlines(y=0, xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.hlines(y=len(table.index), xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.set_xlim(0,len(table.columns))
g.ax_heatmap.set_ylim(0,len(table.index))
g.ax_heatmap.set_xlabel('Gio', fontsize=80)
g.ax_heatmap.set_ylabel('Gprotein Position', fontsize=60)
g.cax.set_position([1.73,0.5, 0.03, 0.5]) 
g.cax.set_title("Conservation",loc="center",pad=30)
for ytic in g.ax_heatmap.get_yticklabels():
    if ytic.get_text() in gp_ann.keys(): # Change color if exist else not
        ytic.set_color(gp_ann[ytic.get_text()])
for xtic in g.ax_heatmap.get_xticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
ann_colors.update(gprot_ann) 
        
handles = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
leg2=g.ax_heatmap.legend(handles, ann_colors,bbox_to_anchor=(1.7,0.06), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=50,title='GPCR_pos/Gprot_pos',title_fontsize=50,columnspacing=0.1,handletextpad =0.1)



g.savefig('./heatmaps/plots/Gio_heatmap_'+fil+'.svg',dpi=600,bbox_inches='tight',format='svg')


###

#Gs
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'},font_scale=5)
g=sns.clustermap(table1,cmap='Reds',row_cluster=False,col_cluster=False,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,vmin=0,vmax=1,linecolor='silver',figsize=(20,20),cbar_kws=dict(use_gridspec=False,pad=10,ticks=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]))
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =50)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 50) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.gs.update(top=1.5,right=1.5)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
ax= divider.append_axes("top", size="12000%", pad=0)
ax2= divider2.append_axes("right", size="12000%", pad=0)
tg.plot.bar(x='BW',y='Gs', color=['r'],ax=ax ,edgecolor='black',width=0.9,legend=False)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_facecolor('white')
ax.set_ylabel('log_odds_ratio',fontsize=50,labelpad=30)
ax.tick_params(axis='y', which='both', labelsize=42)
ax.set_yticks(np.arange(-1,1.2,0.2))
ax.set_ylim(-1.2,1)
tf.plot.barh(x='Gprot_pos',y='Gs',ax=ax2, color=['r'] ,edgecolor='black',fontsize=15,width=0.9,legend=False)
ax2.set_yticklabels([])
ax2.set_ylabel('')
ax2.set_facecolor('white')
ax2.set_xlabel('log_odds_ratio',fontsize=50,labelpad=20)
ax2.tick_params(axis='x', which='both', labelsize=42,rotation=90)
ax2.set_xlim(-1.2,1)
ax2.set_xticks(np.arange(-1,1.2,0.2))
ax2.invert_xaxis()
g.ax_heatmap.hlines(y=0, xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.hlines(y=len(table.index), xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.set_xlim(0,len(table.columns))
g.ax_heatmap.set_ylim(0,len(table.index))
g.ax_heatmap.set_xlabel('Gs', fontsize=80)
g.ax_heatmap.set_ylabel('Gprotein Position', fontsize=60)
g.cax.set_position([1.73,0.5, 0.03, 0.5]) 
g.cax.set_title("Conservation",loc="center",pad=30)
for ytic in g.ax_heatmap.get_yticklabels():
    if ytic.get_text() in gp_ann.keys(): # Change color if exist else not
        ytic.set_color(gp_ann[ytic.get_text()])
for xtic in g.ax_heatmap.get_xticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
handles = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
leg2=g.ax_heatmap.legend(handles, ann_colors,bbox_to_anchor=(1.7,0.06), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=50,title='GPCR_pos/Gprot_pos',title_fontsize=50,columnspacing=0.1,handletextpad =0.1)
        
g.savefig('./heatmaps/plots/Gs_heatmap_'+fil+'.svg',dpi=600,bbox_inches='tight',format='svg')



####Gq11
sns.set(style="ticks",rc={'axes.facecolor':'lightgray'},font_scale=5)
g=sns.clustermap(table3,cmap='Greens',row_cluster=False,col_cluster=False,square=True,xticklabels=1,yticklabels=1,linewidths=0.5,vmin=0,vmax=1,linecolor='silver',figsize=(20,20), cbar_kws={"ticks":[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]})
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =50)
g.ax_heatmap.tick_params(axis='x', which='major', pad=1)
g.ax_heatmap.tick_params(axis='y', which='major', pad=1)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 50) 
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
g.gs.update(top=1.5,right=1.5)
divider = make_axes_locatable(g.ax_col_dendrogram)
divider2 = make_axes_locatable(g.ax_row_dendrogram)
ax= divider.append_axes("top", size="12000%", pad=0)
ax2= divider2.append_axes("right", size="12000%", pad=0)
tg.plot.bar(x='BW',y='Gq11', color=['#008000'],ax=ax ,edgecolor='black',width=0.9,legend=False)
ax.set_xticklabels([])
ax.set_xlabel('')
ax.set_facecolor('white')
ax.set_ylabel('log_odds_ratio',fontsize=50,labelpad=30)
ax.tick_params(axis='y', which='both', labelsize=42)
ax.set_yticks(np.arange(-1,1.2,0.2))
ax.set_ylim(-1.2,1)
tf.plot.barh(x='Gprot_pos',y='Gq11',ax=ax2, color=['#008000'] ,edgecolor='black',fontsize=15,width=0.9,legend=False)
ax2.set_yticklabels([])
ax2.set_ylabel('')
ax2.set_facecolor('white')
ax2.set_xlabel('log_odds_ratio',fontsize=50,labelpad=20)
ax2.tick_params(axis='x', which='both', labelsize=42,rotation=90)
ax2.set_xlim(-1.2,1)
ax2.set_xticks(np.arange(-1,1.2,0.2))
ax2.invert_xaxis()
g.ax_heatmap.hlines(y=0, xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.hlines(y=len(table.index), xmin=0, xmax=len(table.columns), linewidth=5, color='black')
g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=5, color='black')
g.ax_heatmap.set_xlim(0,len(table.columns))
g.ax_heatmap.set_ylim(0,len(table.index))
g.ax_heatmap.set_xlabel('Gq11', fontsize=80)
g.ax_heatmap.set_ylabel('Gprotein Position', fontsize=60)
g.cax.set_position([1.73,0.5, 0.03, 0.5]) 
g.cax.set_title("Conservation",loc="center",pad=30)\

for ytic in g.ax_heatmap.get_yticklabels():
    if ytic.get_text() in gp_ann.keys(): # Change color if exist else not
        ytic.set_color(gp_ann[ytic.get_text()])
for xtic in g.ax_heatmap.get_xticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
handles = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
leg2=g.ax_heatmap.legend(handles, ann_colors,bbox_to_anchor=(1.7,0.06), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=50,title='GPCR_pos/Gprot_pos',title_fontsize=50,columnspacing=0.1,handletextpad =0.1)

g.savefig('./heatmaps/plots/Gq11_heatmap_'+fil+'.svg',dpi=600,bbox_inches='tight',format='svg')







