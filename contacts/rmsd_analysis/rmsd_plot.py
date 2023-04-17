from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from statannot import add_stat_annotation
import sys 
import numpy as np
import pickle

####Plot Heatmaps of rmsd calculation and the violin plots of the distribution of the rmsd.

fs='\t'
fil=sys.argv[1]
#Gs vs Gi pairs

path='./plots/'
df=pd.read_csv(path+fil+".txt",header=None, comment="#", sep="\t")
df.columns=['PDB_ID','PDB_ID1','rmsd']

df = df.replace(r'\n',' ', regex=True) 
df['PDB_ID'] = df['PDB_ID'].astype(str).str.strip()
df['PDB_ID1'] = df['PDB_ID1'].astype(str).str.strip()

df1=pd.read_csv("./use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
ds=df1[['PDB_ID','Receptor_Gene_name', 'Gprotein_Gene_name','Antibody','Class']].drop_duplicates()
ds['name']=ds['Receptor_Gene_name']+"_"+ ds['Gprotein_Gene_name']+"_"+ds['PDB_ID']
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
ds["gfam"] = ds["Gprotein_Gene_name"].map(equiv)
dp=ds[['PDB_ID','name']].drop_duplicates()
dp['PDB_ID'] = dp['PDB_ID'].astype(str)
dp['name'] = dp['name'].astype(str)
dp=dp.reset_index(drop=True)


mydict = dict(zip(dp.PDB_ID, dp.name))
df2=df.replace({"PDB_ID": mydict})
df2=df2.replace({"PDB_ID1": mydict})
df2[['GPCR1', 'Gprotein1', 'pdb1']] = df2['PDB_ID'].str.split('_', expand=True)
df2[['GPCR2', 'Gprotein2', 'pdb2']] = df2['PDB_ID1'].str.split('_', expand=True)

equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df2["gfam1"] = df2["Gprotein1"].map(equiv)
df2["gfam2"] = df2["Gprotein2"].map(equiv)
df2=df2.dropna()
df2.to_csv(path+fil+".tsv",index=None,header=True,sep='\t')

#draw
ds2=df2[['PDB_ID', 'gfam1']].drop_duplicates()
ds3=df2[['PDB_ID1', 'gfam2']].drop_duplicates()
ds2.columns=['PDB_ID','gfam']
ds3.columns=['PDB_ID','gfam']
ds3=pd.concat([ds2,ds3]).drop_duplicates()
ds3.set_index('PDB_ID',inplace=True)
ds3 = ds3.pop("gfam")
ds3=ds3.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}
rs = ds3.map(lut)
#class
dc=df2[['PDB_ID', 'GPCR1']].drop_duplicates()
dc1=df2[['PDB_ID1', 'GPCR2']].drop_duplicates()
dc.columns=['PDB_ID','GPCR']
dc1.columns=['PDB_ID','GPCR']
dc2=pd.concat([dc,dc1]).drop_duplicates()

####colors
infile = open('./use_file/GPCR_class.pickle', 'rb')
clas = pickle.load(infile)

# create the customized color map gfamily
ds1=ds[['name', 'gfam']].drop_duplicates()
ds1.set_index('name',inplace=True)
ds1 = ds1.pop("gfam")
ds1=ds1.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}
rs = ds1.map(lut)
# create the customized color map gpcr class

clas1=ds[['name', 'Class']].drop_duplicates().set_index('name')
clas1=clas1.astype('object')
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}
rs1 = clas1['Class'].map(lut1)

infile3 = open('./use_file/Coupling_assay.pickle', 'rb')
coup = pickle.load(infile3)
rs3 = pd.concat([rs,rs1],axis=1)
co=[]
for index, row in rs3.iterrows():
    key = index.split('_')[0]  
    
    # check if the key exists in coup
    if key in coup:
        v = coup[key]
        co.append((index, row['Class'], row['gfam'], v[0], v[1], v[2], v[3]))
    else:
        # add a tuple with None values for the missing data
        co.append((index, row['Class'], row['gfam'], None, None, None, None))
coup1=pd.DataFrame(co,columns=['GPCR','Class','gfam','Gs','Gio','Gq11','G1213']).set_index('GPCR')
rs3 = coup1
lut.update(lut1)
#coup={'G1213':'#800080','Only GtoPdb':'#cfccc6','Only TGF or GEMTA':'#FF8C00','Not coupled':'#808080'}
cou={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}
lut.update(cou)

ab=ds[['name', 'Antibody']].drop_duplicates().set_index('name')
a1={'No':'#f2f2f2'}
ab["antibody"] = ab["Antibody"].map(a1)
ab=ab.fillna('#000000')

ab1={'Antibody':'#000000','No antibody':'#f2f2f2'}
lut.update(ab1) 
ab=ab.drop(['Antibody'], axis=1)
rs3 = pd.concat([rs3,ab],axis=1)
rs3 = rs3[ ['gfam'] + [ col for col in rs3.columns if col != 'gfam' ] ]


#classA
clas1=clas1.reset_index()
mydict1 = dict(zip(df1['PDB_ID'],clas1['Class']))
df2["class1"] = df2["pdb1"].map(mydict1)
df2["class2"] = df2["pdb2"].map(mydict1)
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
#all class

def plot_clustermap(table,fil=str(),name=str()):    
    mask=0
    sns.set(font_scale=6.5)
    num_rows, num_cols = table.shape
    figsize = (num_cols //4, num_rows //4)
    sns.set(style="ticks",rc={'axes.facecolor':'lightgray'})
    g=sns.clustermap(table,method='ward',mask=table==mask,col_colors=rs3,row_colors=rs3,tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.08, 0.08),cmap='Reds',square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver',colors_ratio=(0.007,0.007),figsize=figsize)
    g.ax_heatmap.set_xlabel('Position', fontsize=60, fontweight='bold')
    g.ax_heatmap.set_ylabel('GPCR', fontsize=60, fontweight='bold')
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize =18, fontweight='bold')
    g.ax_heatmap.tick_params(axis='x', which='major', pad=2)
    g.ax_heatmap.tick_params(axis='y', which='major', pad=2)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = 18, fontweight='bold') 
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    #g.gs.update(top=0.1,right=0.05)
    g.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table.columns), linewidth=4, color='black')
    g.ax_heatmap.hlines(y=len(table.index), xmin=-0.5, xmax=len(table.columns), linewidth=5, color='black')
    g.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=5, color='black')
    g.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=5, color='black')
    g.cax.set_position([1.05, 0.45, 0.01, 0.2]) 
    #g.cax.set_yticklabels(g.cax.get_ymajorticklabels(), fontsize =22)
    g.cax.tick_params(labelsize=60)
    g.ax_row_colors.tick_params(labelsize=35)
    g.ax_col_colors.tick_params(labelsize=35)
    g.cax.set_title("RMSD",loc="center",pad=20,fontsize=70)
    #gfamily and class
    handles = [Patch(facecolor=lut[name]) for name in lut]
    leg=g.ax_heatmap.legend(handles, lut,bbox_to_anchor=(0.6,-0.028), handleheight=4, handlelength=4,title_fontsize=40, bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=5,title='Gprotein family/GPCR class', prop={'weight': 'bold'})
    for label in leg.get_texts():
        label.set_fontsize(35)
    g.savefig('./plots/'+fil+'_htmp_'+name+'.svg',bbox_inches='tight',dpi=300)

plot_clustermap(table1,fil,'all')
plot_clustermap(table2,fil,'classA') 


#violin plots full
dg=df2[(df2['gfam1']=='Gs') & (df2['gfam2']=='Gs')]
di=df2[(df2['gfam1']=='Gio') & (df2['gfam2']=='Gio')]
dq=df2[(df2['gfam1']=='Gq11') & (df2['gfam2']=='Gq11')]
dgi=pd.concat([dg,di])
dgi.to_csv("./plots/violin_Gs_Gi_all_"+fil+".tsv",index=None,header=True,sep='\t')


#class A
#violin plots 
dgA=df3[(df3['gfam1']=='Gs') & (df3['gfam2']=='Gs')]
diA=df3[(df3['gfam1']=='Gio') & (df3['gfam2']=='Gio')]
dqA=df3[(df3['gfam1']=='Gq11') & (df3['gfam2']=='Gq11')]
dgiA=pd.concat([dgA,diA])
dgiA.to_csv("./plots/violin_Gs_Gi_classA_"+fil+".tsv",index=None,header=True,sep='\t')

def plot_violinplot(table,t1,t2,fil=str(),name=str()):    
    sns.set(style="ticks",rc={'axes.facecolor':'white'},font_scale=2)
    fig , ax = plt.subplots(figsize=(10,10))
    ax = sns.violinplot(x=table["gfam1"],y=table["rmsd"],palette=lut, dodge=False)
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    pairs=[('Gs', 'Gio')]
    stat,pval=stats.ranksums(t1["rmsd"], t2["rmsd"])
    pvalues = [pval]
    test_short_name = 'Ranksum'
    ax.set_ylim([table['rmsd'].min()-2, table['rmsd'].max()+5])  
    ax.set_xlabel("G-protein", fontsize=50)
    ax.set_ylabel("RMSD",fontsize=50)
    ax.tick_params(labelsize=50)
    t= add_stat_annotation(ax, data=table, x='gfam1', y='rmsd',test_short_name=test_short_name,
                                   box_pairs=pairs, perform_stat_test=False,
                                   test=None, text_format='full',pvalues=pvalues,
                                   loc='inside', verbose=2)
    plt.savefig('./plots/'+fil+'_vlp_'+name+'.svg',bbox_inches='tight')
    with open('./plots/'+fil+'_stats_'+name+'.txt', 'w') as f:
        print('Gprotein',fs,'Median',fs,'std',fs,'sem',file=f)
        print ('Gs',fs,np.median(t1["rmsd"]),fs, np.std(t1["rmsd"]),fs, stats.sem(t1["rmsd"]),file=f)
        print ('Gio',fs,np.median(t2["rmsd"]), np.std(t2["rmsd"]), stats.sem(t2["rmsd"]),file=f)
        #print ('Gq11',fs,np.median(dq["rmsd"]), np.std(dq["rmsd"]), stats.sem(dq["rmsd"]),file=f)

plot_violinplot(dgi,dg,di,fil,'all')
plot_violinplot(dgiA,dgA,diA,fil,'classA')

