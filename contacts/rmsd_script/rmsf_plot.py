#!/data/SW/anaconda3/envs/myenv/bin/python

from matplotlib.patches import Patch
from scipy import stats
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#files
df1=pd.read_csv("gprot_all_gpcr_all_rmsf.tsv", comment="#", sep="\t").dropna(axis=1,how='all')
df2=pd.read_csv("gprot_all_gpcr_all12_rmsf.tsv", comment="#", sep="\t").dropna(axis=1,how='all')
df3=pd.read_csv("gprot_cons_gpcr_all_rmsf.tsv", comment="#", sep="\t").dropna(axis=1,how='all')
df4=pd.read_csv("gprot_cons_gpcr_all12_rmsf.tsv", comment="#", sep="\t").dropna(axis=1,how='all')

#positions
d1=pd.read_csv("gprot_pos_set_all.tsv", sep="\t").dropna()
d1[['G', 'Gprot','num']] = d1['CGN'].str.split('.', 2, expand=True)
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
d1.Gprot = d1.Gprot.astype("category")
d1.Gprot.cat.set_categories(sorter1, inplace=True)
d1.sort_values(["Gprot","CGN"],inplace=True)
d1=d1.reset_index()
d1=d1.reset_index()
#patches
#gprot ann
gprot_ann={'Helix':'deepskyblue','Sheet':'crimson','Loop':'darkorchid'}
gprotp=d1['Gprot'].to_frame().drop_duplicates()
gprotp1 = d1[~d1['Gprot'].str.isupper()]
gprotp1['col']='darkorchid'
gprotp2 = gprotp[gprotp['Gprot'].str.isupper() & gprotp['Gprot'].str.startswith('H')]
gprotp3 = gprotp[gprotp['Gprot'].str.isupper() & gprotp['Gprot'].str.startswith('S')]
gprotp2['col']='deepskyblue'
gprotp3['col']='crimson'
gp=pd.concat([gprotp1,gprotp2,gprotp3])
gp_ann = dict(zip(gp.Gprot, gp.col))
min1=d1.groupby(['Gprot'])['level_0'].min().reset_index().dropna()
max1=d1.groupby(['Gprot'])['level_0'].max().reset_index().dropna()
pos = min1.merge(max1,left_on=['Gprot'],right_on=['Gprot'],how="left")
posd=pos.set_index('Gprot').transpose().to_dict(orient='list')

        
#mean and std and plot gprot_all_gpcr_all

dg=df1[(df1['gfam1']=='Gs')]
dg=dg[dg.isin(['6p9x']).any(axis=1)]
di=df1[(df1['gfam1']=='Gio')]
di=di[di.isin(['7y15']).any(axis=1)]

gsmean=dg.mean().reset_index()
gsstd=dg.std().reset_index()
giomean=di.mean().reset_index()
giostd=di.std().reset_index()
gios = giomean.merge(giostd,left_on=['index'],right_on=['index'],how="left")
ggs = gsmean.merge(gsstd,left_on=['index'],right_on=['index'],how="left")
stat,pval=stats.wilcoxon(ggs["0_x"], gios["0_x"])
#plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})
plt.figure(figsize=(15,8))
plt.rcParams['figure.facecolor'] = 'white'
plt.plot(ggs['index'], ggs['0_x'], 'r-', label='Gs')
plt.fill_between(ggs['index'], ggs['0_x'] - ggs['0_y'], ggs['0_x'] + ggs['0_y'], color='r', alpha=0.2)
plt.plot(gios['index'], gios['0_x'], 'b-', label='Gio')
plt.fill_between(gios['index'], gios['0_x'] - gios['0_y'], gios['0_x'] + gios['0_y'], color='b', alpha=0.2)
plt.hlines(y=-2.1, xmin=0, xmax=81, colors='dimgrey', linestyles='-', lw=4)
for key in gp_ann.keys():
    for key1 in posd.keys():
        if key==key1:
            plt.hlines(-2, posd[key1][0]-0.25,posd[key1][1]+0.5, colors=gp_ann[key], linestyles='solid', linewidths=30)
            plt.text((posd[key1][0]-3+posd[key1][1])/2, -3.85,key1, size=20, color='black')
plt.xlim(-1,81)
plt.xticks([])
plt.yticks(np.arange(-1,17, 1),fontsize=20)
plt.ylim(-3,16)
plt.xlabel('Secondary Structure',fontsize=20,labelpad=20)
plt.ylabel('Mean RMSF±STD',fontsize=20)
gis={'Gs':'r','Gio':'b'}
gprot_ann.update(gis) 
plt.text(55,15,'Wilcoxon pvalue-'+str('{:.2e}'.format(pval)), size=20, color='black')
handles = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
plt.legend(handles, gprot_ann,bbox_to_anchor=(0.73,0.055), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=20,title='GPCR_pos/Gprot_pos',title_fontsize=20,columnspacing=0.1,handletextpad =0.1)
plt.savefig("./plots/rmsf_gprot_all_gpcr_all.svg",bbox_inches='tight', dpi=600)

   
#mean and std and plot gprot_all_gpcr_all12
dg=df2[(df2['gfam1']=='Gs')]
dg=dg[dg.isin(['6p9x']).any(axis=1)]
di=df2[(df2['gfam1']=='Gio')]
di=di[di.isin(['7y15']).any(axis=1)]
gsmean=dg.mean().reset_index()
gsstd=dg.std().reset_index()
giomean=di.mean().reset_index()
giostd=di.std().reset_index()
gios = giomean.merge(giostd,left_on=['index'],right_on=['index'],how="left")
ggs = gsmean.merge(gsstd,left_on=['index'],right_on=['index'],how="left")
stat1,pval1=stats.wilcoxon(ggs["0_x"], gios["0_x"])
#plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})
plt.figure(figsize=(15,8))
plt.rcParams['figure.facecolor'] = 'white'
plt.plot(ggs['index'], ggs['0_x'], 'r-', label='Gs')
plt.fill_between(ggs['index'], ggs['0_x'] - ggs['0_y'], ggs['0_x'] + ggs['0_y'], color='r', alpha=0.2)
plt.plot(gios['index'], gios['0_x'], 'b-', label='Gio')
plt.fill_between(gios['index'], gios['0_x'] - gios['0_y'], gios['0_x'] + gios['0_y'], color='b', alpha=0.2)
plt.hlines(y=-2.1, xmin=0, xmax=81, colors='dimgrey', linestyles='-', lw=4)
for key in gp_ann.keys():
    for key1 in posd.keys():
        if key==key1:
            plt.hlines(-2, posd[key1][0]-0.25,posd[key1][1]+0.5, colors=gp_ann[key], linestyles='solid', linewidths=30)
            plt.text((posd[key1][0]-3+posd[key1][1])/2, -3.85,key1, size=20, color='black')
plt.xlim(-1,81)
plt.xticks([])
plt.yticks(np.arange(-1,17, 1),fontsize=20)
plt.ylim(-3,16)
plt.xlabel('Secondary Structure',fontsize=20,labelpad=20)
plt.ylabel('Mean RMSF±STD',fontsize=20)
gis={'Gs':'r','Gio':'b'}
gprot_ann.update(gis) 
plt.text(55,15,'Wilcoxon pvalue-'+str('{:.2e}'.format(pval1)), size=20, color='black')
handles = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
plt.legend(handles, gprot_ann,bbox_to_anchor=(0.73,0.055), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=20,title='GPCR_pos/Gprot_pos',title_fontsize=20,columnspacing=0.1,handletextpad =0.1)
plt.savefig("./plots/rmsf_gprot_all_gpcr_all12.svg",bbox_inches='tight', dpi=600)

#positions
d1=pd.read_csv("gprot_pos_set_cons.tsv", sep="\t").dropna()
d1[['G', 'Gprot','num']] = d1['CGN'].str.split('.', 2, expand=True)
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
d1.Gprot = d1.Gprot.astype("category")
d1.Gprot.cat.set_categories(sorter1, inplace=True)
d1.sort_values(["Gprot","CGN"],inplace=True)
d1=d1.reset_index()
d1=d1.reset_index()
#patches
#gprot ann
gprot_ann={'Helix':'deepskyblue','Sheet':'crimson','Loop':'darkorchid'}
gprotp=d1['Gprot'].to_frame().drop_duplicates()
gprotp1 = d1[~d1['Gprot'].str.isupper()]
gprotp1['col']='darkorchid'
gprotp2 = gprotp[gprotp['Gprot'].str.isupper() & gprotp['Gprot'].str.startswith('H')]
gprotp3 = gprotp[gprotp['Gprot'].str.isupper() & gprotp['Gprot'].str.startswith('S')]
gprotp2['col']='deepskyblue'
gprotp3['col']='crimson'
gp=pd.concat([gprotp1,gprotp2,gprotp3])
gp_ann = dict(zip(gp.Gprot, gp.col))
min1=d1.groupby(['Gprot'])['level_0'].min().reset_index().dropna()
max1=d1.groupby(['Gprot'])['level_0'].max().reset_index().dropna()
pos = min1.merge(max1,left_on=['Gprot'],right_on=['Gprot'],how="left")
posd=pos.set_index('Gprot').transpose().to_dict(orient='list')





       
#mean and std and plot gprot_cons_gpcr_all
dg=df3[(df3['gfam1']=='Gs')]
dg=dg[dg.isin(['6p9x']).any(axis=1)]
di=df3[(df3['gfam1']=='Gio')]
di=di[di.isin(['7y15']).any(axis=1)]
gsmean=dg.mean().reset_index()
gsstd=dg.std().reset_index()
giomean=di.mean().reset_index()
giostd=di.std().reset_index()
gios = giomean.merge(giostd,left_on=['index'],right_on=['index'],how="left")
ggs = gsmean.merge(gsstd,left_on=['index'],right_on=['index'],how="left")
stat,pval=stats.wilcoxon(ggs["0_x"], gios["0_x"])
#plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})
plt.figure(figsize=(15,8))
plt.rcParams['figure.facecolor'] = 'white'
plt.plot(ggs['index'], ggs['0_x'], 'r-', label='Gs')
plt.fill_between(ggs['index'], ggs['0_x'] - ggs['0_y'], ggs['0_x'] + ggs['0_y'], color='r', alpha=0.2)
plt.plot(gios['index'], gios['0_x'], 'b-', label='Gio')
plt.fill_between(gios['index'], gios['0_x'] - gios['0_y'], gios['0_x'] + gios['0_y'], color='b', alpha=0.2)
plt.hlines(y=-2.1, xmin=0, xmax=63, colors='dimgrey', linestyles='-', lw=4)
for key in gp_ann.keys():
    for key1 in posd.keys():
        if key==key1:
            plt.hlines(-2, posd[key1][0]-0.25,posd[key1][1]+0.5, colors=gp_ann[key], linestyles='solid', linewidths=30)
            plt.text((posd[key1][0]-3+posd[key1][1])/2, -3.85,key1, size=20, color='black')
plt.xlim(-1,63)
plt.xticks([])
plt.yticks(np.arange(-1,17, 1),fontsize=20)
plt.ylim(-3,16)
plt.xlabel('Secondary Structure',fontsize=20,labelpad=20)
plt.ylabel('Mean RMSF±STD',fontsize=20)
gis={'Gs':'r','Gio':'b'}
gprot_ann.update(gis) 
plt.text(35,15,'Wilcoxon pvalue-'+str('{:.2e}'.format(pval)), size=20, color='black')
handles = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
plt.legend(handles, gprot_ann,bbox_to_anchor=(0.73,0.055), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=20,title='GPCR_pos/Gprot_pos',title_fontsize=20,columnspacing=0.1,handletextpad =0.1)
plt.savefig("./plots/rmsf_gprot_cons_gpcr_all.svg",bbox_inches='tight', dpi=600)

   
#mean and std and plot gprot_all_gpcr_all12
dg=df4[(df4['gfam1']=='Gs')]
di=df4[(df4['gfam1']=='Gio')]
dg=dg[dg.isin(['6p9x']).any(axis=1)]
di=di[di.isin(['7y15']).any(axis=1)]
gsmean=dg.mean().reset_index()
gsstd=dg.std().reset_index()
giomean=di.mean().reset_index()
giostd=di.std().reset_index()
gios = giomean.merge(giostd,left_on=['index'],right_on=['index'],how="left")
ggs = gsmean.merge(gsstd,left_on=['index'],right_on=['index'],how="left")
stat1,pval1=stats.wilcoxon(ggs["0_x"], gios["0_x"])
#plot
sns.set(style="ticks",rc={'axes.facecolor':'white'})
plt.figure(figsize=(15,8))
plt.rcParams['figure.facecolor'] = 'white'
plt.plot(ggs['index'], ggs['0_x'], 'r-', label='Gs')
plt.fill_between(ggs['index'], ggs['0_x'] - ggs['0_y'], ggs['0_x'] + ggs['0_y'], color='r', alpha=0.2)
plt.plot(gios['index'], gios['0_x'], 'b-', label='Gio')
plt.fill_between(gios['index'], gios['0_x'] - gios['0_y'], gios['0_x'] + gios['0_y'], color='b', alpha=0.2)
plt.hlines(y=-2.1, xmin=0, xmax=63, colors='dimgrey', linestyles='-', lw=4)
for key in gp_ann.keys():
    for key1 in posd.keys():
        if key==key1:
            plt.hlines(-2, posd[key1][0]-0.25,posd[key1][1]+0.5, colors=gp_ann[key], linestyles='solid', linewidths=30)
            plt.text((posd[key1][0]-3+posd[key1][1])/2, -3.85,key1, size=20, color='black')
plt.xlim(-1,63)
plt.xticks([])
plt.yticks(np.arange(-1,17, 1),fontsize=20)
plt.ylim(-3,16)
plt.xlabel('Secondary Structure',fontsize=20,labelpad=20)
plt.ylabel('Mean RMSF±STD',fontsize=20)
gis={'Gs':'r','Gio':'b'}
gprot_ann.update(gis) 
plt.text(35,15,'Wilcoxon pvalue-'+str('{:.2e}'.format(pval1)), size=20, color='black')
handles = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
plt.legend(handles, gprot_ann,bbox_to_anchor=(0.73,0.055), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=20,title='GPCR_pos/Gprot_pos',title_fontsize=20,columnspacing=0.1,handletextpad =0.1)
plt.savefig("./plots/rmsf_gprot_cons_gpcr_all12.svg",bbox_inches='tight', dpi=600)


