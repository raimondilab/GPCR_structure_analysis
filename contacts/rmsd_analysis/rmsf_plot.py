from matplotlib.patches import Patch
from scipy import stats
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sys
#################### Plotting RMSF distribution of Gs vs Gio
nom=sys.argv[1]
df=pd.read_csv("./plots/violin_Gs_Gi_all_bioalign_"+nom+".tsv", comment="#", sep="\t")
#violin plots full
dg=df[(df['gfam1']=='Gs') & (df['gfam2']=='Gs')]
di=df[(df['gfam1']=='Gio') & (df['gfam2']=='Gio')]
#representative coupling lowest rmsd mean
dg=((dg.groupby(['GPCR1'])['rmsd'].mean()+dg.groupby(['GPCR2'])['rmsd'].mean())/2).to_frame()
di=((di.groupby(['GPCR1'])['rmsd'].mean()+di.groupby(['GPCR2'])['rmsd'].mean())/2).to_frame()
low_gs=dg['rmsd'].idxmin()
low_gi=di['rmsd'].idxmin()
with open('./plots/'+nom+'_rmsf_'+'lowest_rmsd.txt', 'w') as j:
    print(nom,'\t',low_gs,'\t',low_gi,file=j)


#files
df1=pd.read_csv("./plots/"+nom+"_rmsf.tsv", comment="#", sep="\t").dropna(axis=1,how='all')
#positions
d1=pd.read_csv("./pos/gprot_pos_set_all.tsv", sep="\t").dropna()
d1[['Gprot','num']] = d1['Gprot_pos'].str.split('.', 2, expand=True)
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
d1.Gprot = d1.Gprot.astype("category")
d1.Gprot.cat.set_categories(sorter1, inplace=True)
d1['num']=d1['num'].astype(int)
d1.sort_values(by=['Gprot',"num"],ascending=True,inplace=True)
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
gs=df1[(df1['gfam1']=='Gs')]
gs=gs[gs.isin([low_gs]).any(axis=1)]
gio=df1[(df1['gfam1']=='Gio')]
gio=gio[gio.isin([low_gi]).any(axis=1)]
gs=gs.drop(['rmsd'], axis=1)
gio=gio.drop(['rmsd'], axis=1)
gsmean=gs.mean().reset_index()
gsstd=gs.std().reset_index()
giomean=gio.mean().reset_index()
giostd=gio.std().reset_index()
gios = giomean.merge(giostd,left_on=['index'],right_on=['index'],how="left")
ggs = gsmean.merge(gsstd,left_on=['index'],right_on=['index'],how="left")
stat,pval=stats.wilcoxon(ggs["0_x"], gios["0_x"])

def plot_rmsf(t1,t2,name=str()):
    sns.set(style="ticks",rc={'axes.facecolor':'white'})
    plt.figure(figsize=(15,8))
    ym=round(t2['0_x'].max()+t2['0_y'].max())+2
    xm=len(t1.index)
    plt.rcParams['figure.facecolor'] = 'white'
    plt.plot(t1['index'], t1['0_x'], 'r-', label='Gs')
    plt.fill_between(t1['index'], t1['0_x'] - t1['0_y'], t1['0_x'] + t1['0_y'], color='r', alpha=0.2)
    plt.plot(t2['index'], t2['0_x'], 'b-', label='Gio')
    plt.fill_between(t2['index'], t2['0_x'] - t2['0_y'], t2['0_x'] + t2['0_y'], color='b', alpha=0.2)
    plt.hlines(y=-2.1, xmin=0, xmax=xm, colors='dimgrey', linestyles='-', lw=4)
    for key in gp_ann.keys():
        for key1 in posd.keys():
            if key==key1:
                plt.hlines(-2, posd[key1][0]-0.25,posd[key1][1]+0.5, colors=gp_ann[key], linestyles='solid', linewidths=30)
                plt.text((posd[key1][0]-3+posd[key1][1])/2, -3.85,key1, size=18, color='black', weight='bold')
                plt.xlim(-1,xm)
                plt.xticks([])
                plt.yticks(np.arange(-1,ym, 1),fontsize=25, weight='bold')
                plt.ylim(-3,ym)
                plt.xlabel('Secondary Structure',fontsize=22,labelpad=20, weight='bold')
                plt.ylabel('Mean RMSF ±STD',fontsize=22, weight='bold')
                gis={'Gs':'r','Gio':'b'}
                gprot_ann.update(gis) 
                plt.text(xm/1.8,ym-1,'Wilcoxon pvalue-'+str('{:.2e}'.format(pval)), weight='bold', size=20, color='black')
                handles = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
                plt.legend(handles, gprot_ann,bbox_to_anchor=(0.73,0.055),bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=10,fontsize=20,title='GPCR_pos/Gprot_pos',title_fontsize=20,columnspacing=0.1,handletextpad =0.1)
                plt.savefig("./plots/"+"rmsf_"+name+".svg",bbox_inches='tight', dpi=600)
plot_rmsf(ggs,gios,nom)
 