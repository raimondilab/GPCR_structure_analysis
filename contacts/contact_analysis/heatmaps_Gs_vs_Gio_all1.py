
from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.preprocessing import MaxAbsScaler
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle

#####Plot heatmaps with different filters of contact pairs or GPCR and Gprotein separatley for all structures





#Gs vs Gi pairs
df=pd.read_csv("./use_file/all_gpcrs.tsv", comment="#", sep="\t")
infile = open('./use_file/GPCR_class.pickle', 'rb')
clas = pickle.load(infile)
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein"].map(equiv)
df['comb']=df["BW"]+"-"+df["Gprot_pos"]
df["Class"] = df["GPCR"].map(clas)
infile1 = open('./use_file/GPCR_family.pickle', 'rb')
fam = pickle.load(infile1)
df["GPCR_family"] = df['GPCR'].map(fam)

###annotation
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
t = df[["BW", "hex"]]
mydict = dict(zip(t.BW, t.hex))




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


###duplicate gpcrs
dupl = df1[['GPCR', 'gfam']].drop_duplicates()
dupl_mask = dupl['GPCR'].duplicated()
gpcr_dupl = dupl[dupl_mask]
for index,row in gpcr_dupl.iterrows():
    df1['GPCR']=np.where(((df1['GPCR']==row['GPCR']) & (df1['gfam']==row['gfam'])),row['GPCR']+"-"+row['gfam'],df1['GPCR'])

###gpcr and gprot sort
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
ds1=df1[['Gprot_struct','Gprot_pos']].drop_duplicates()
ds1[['G','num']] = df1['Gprot_pos'].str.split('.', expand=True)
ds1['num']=ds1['num'].astype(int)
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2=ds1["Gprot_pos"].drop_duplicates().tolist()
#gprot
table1 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1 = table1.reindex(columns=lr2)
#gpcr
table12 = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12= table12.reindex(columns=lr1)
mask =0
#02
df1gc.Structural = df1gc.Structural.astype("category")
df1gc.Structural.cat.set_categories(sorter, inplace=True)
df1gc.sort_values(["Structural","BW"],inplace=True)
lr1g=df1gc["BW"].drop_duplicates().tolist()
df1gp.Gprot_struct = df1gp.Gprot_struct.astype("category")
df1gp.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1=df1gp[['Gprot_struct','Gprot_pos']].drop_duplicates()
ds1[['G','num']] = df1['Gprot_pos'].str.split('.', expand=True)
ds1['num']=ds1['num'].astype(int)
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2g=ds1["Gprot_pos"].drop_duplicates().tolist()
#gprot
table1g = pd.pivot_table(df1gp ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
table1g = table1g.reindex(columns=lr2g)
#gpcr
table12g = pd.pivot_table(df1gc ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12g= table12g.reindex(columns=lr1g)
mask =0

######pair
gpcrsp=df1.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrsp = gpcrsp[(gpcrsp['Total_count']/len(set(df1['GPCR'])) > 0.2)]
df1gcp=gpcrsp.merge(df1,left_on=['comb'],right_on=['comb'],how="left")


df1.sort_values(["Structural","comb"],inplace=True)
lr1p=df1["comb"].drop_duplicates().tolist()

#comb
table1p = pd.pivot_table(df1 ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
#table1= table1.reindex(lr2)
table1p= table1p.reindex(columns=lr1p)
#02
df1gcp.Structural = df1gcp.Structural.astype("category")
df1gcp.Structural.cat.set_categories(sorter, inplace=True)
df1gcp.sort_values(["Structural","comb"],inplace=True)
lr1gp=df1gcp["comb"].drop_duplicates().tolist()
#gpcrp
table12gp = pd.pivot_table(df1gcp ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
table12gp= table12gp.reindex(columns=lr1gp)
mask =0



####classA
#gprot
df1A=df1.loc[df1['Class'] == 'classA']
gpcrsA=df1A.groupby(['BW'])['GPCR'].nunique().reset_index(name='Total_count')
gprotsA=df1A.groupby(['Gprot_pos'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrsA = gpcrsA[(gpcrsA['Total_count']/len(set(df1A['GPCR'])) > 0.2)]
gprotsA = gprotsA[(gprotsA['Total_count']/len(set(df1A['GPCR'])) > 0.2)]
df1gcA=gpcrsA.merge(df1A,left_on=['BW'],right_on=['BW'],how="left")
df1gpA=gprotsA.merge(df1A,left_on=['Gprot_pos'],right_on=['Gprot_pos'],how="left")
df1A.Structural = df1A.Structural.astype("category")
df1A.Structural.cat.set_categories(sorter, inplace=True)
df1A.sort_values(["Structural","BW"],inplace=True)
lr1A=df1A["BW"].drop_duplicates().tolist()
df1A.Gprot_struct = df1A.Gprot_struct.astype("category")
df1A.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1A=df1A[['Gprot_struct','Gprot_pos']].drop_duplicates()
ds1A[['G','num']] = df1A['Gprot_pos'].str.split('.', expand=True)
ds1A['num']=ds1A['num'].astype(int)
ds1A.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2A=ds1A["Gprot_pos"].drop_duplicates().tolist()
#gprot
table1A = pd.pivot_table(df1A ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1A = table1A.reindex(columns=lr2A)
table1A = pd.pivot_table(df1A ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
#table1= table1.reindex(lr2)
table1A = table1A.reindex(columns=lr2A)
table1A=table1A.dropna(axis=1, how='all')
#gpcr
table12A = pd.pivot_table(df1A ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12A= table12A.reindex(columns=lr1A)
table12A=table12A.dropna(axis=1, how='all')
mask =0
table12A=table12A.drop(['ICL1','Cterm','H8','ICL2','ICL3'], axis=1)
lr1A.remove('ICL1')
lr1A.remove('ICL2')
lr1A.remove('ICL3')
lr1A.remove('H8')
lr1A.remove('Cterm')
#gprot
###02
df1gcA.Structural = df1gcA.Structural.astype("category")
df1gcA.Structural.cat.set_categories(sorter, inplace=True)
df1gcA.sort_values(["Structural","BW"],inplace=True)
lr1gA=df1gcA["BW"].drop_duplicates().tolist()
df1gpA.Gprot_struct = df1gpA.Gprot_struct.astype("category")
df1gpA.Gprot_struct.cat.set_categories(sorter1, inplace=True)
ds1A=df1gpA[['Gprot_struct','Gprot_pos']].drop_duplicates()
ds1A[['G','num']] = df1gpA['Gprot_pos'].str.split('.', expand=True)
ds1A['num']=ds1A['num'].astype(int)
ds1A.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
lr2gA=ds1A["Gprot_pos"].drop_duplicates().tolist()
table1gA = pd.pivot_table(df1gpA ,values='val', index=['GPCR'],columns=['Gprot_pos']).fillna(0)
table1gA = table1gA.reindex(columns=lr2g)
table1gA=table1gA.dropna(axis=1, how='all')
##02
table12gA = pd.pivot_table(df1gcA ,values='val', index=['GPCR'],columns=['BW']).fillna(0)
table12gA= table12gA.reindex(columns=lr1g)
table12gA=table12gA.dropna(axis=1, how='all')
table12gA=table12gA.drop(['ICL3'], axis=1)
lr1gA.remove('ICL3')

##############pair classA

df1Ap=df1A[(df1A.BW != 'ICL1') & (df1A.BW != 'Cterm') & (df1A.BW != 'ICL2') & (df1A.BW != 'ICL3') & (df1A.BW != 'H8')]
gpcrsAp=df1Ap.groupby(['comb'])['GPCR'].nunique().reset_index(name='Total_count')
gpcrsAp = gpcrsAp[(gpcrsAp['Total_count']/len(set(df1Ap['GPCR'])) > 0.2)]
df1gcAp=gpcrsAp.merge(df1Ap,left_on=['comb'],right_on=['comb'],how="left")
df1Ap.Structural = df1Ap.Structural.astype("category")
df1Ap.Structural.cat.set_categories(sorter, inplace=True)
df1Ap.sort_values(["Structural","comb"],inplace=True)
lr1Ap=df1Ap["comb"].drop_duplicates().tolist()

table1Ap = pd.pivot_table(df1Ap ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
table1Ap = table1Ap.reindex(columns=lr1Ap)
table1Ap=table1Ap.dropna(axis=1, how='all')
df1gcAp.Structural = df1gcAp.Structural.astype("category")
df1gcAp.Structural.cat.set_categories(sorter, inplace=True)
df1gcAp.sort_values(["Structural","comb"],inplace=True)
lr1gAp=df1gcAp["comb"].drop_duplicates().tolist()
table1gAp = pd.pivot_table(df1gcAp ,values='val', index=['GPCR'],columns=['comb']).fillna(0)
table1gAp = table1gAp.reindex(columns=lr1gAp)
table1gAp=table1gAp.dropna(axis=1, how='all')


# create the customized color map gfamily
ds=df1[['GPCR', 'gfam']].drop_duplicates()
ds.set_index('GPCR',inplace=True)
ds = ds.pop("gfam")
ds=ds.astype('object')
lut = {'Gio':'b','Gq11':'#008000','Gs':'r','G1213':'#800080'}
rs = ds.map(lut)
# create the customized color map gpcr class

clas1=df1[['GPCR', 'Class']].drop_duplicates().set_index('GPCR')
clas1=clas1.astype('object')
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}
rs1 = clas1['Class'].map(lut1)
# create the customized color map gpcr family
fam1=df1[['GPCR', 'GPCR_family']].drop_duplicates().set_index('GPCR')
fam1['GPCR_family'] = fam1['GPCR_family'].str.replace(r'receptors', '')
fam1['GPCR_family'] = fam1['GPCR_family'].str.replace(r'receptor', '')

fam1=fam1.astype('object')
lut2 = dict(zip(set(fam1['GPCR_family']), sns.husl_palette(len(set(fam1['GPCR_family'])),h=.5)))
rs2= fam1['GPCR_family'].map(lut2)
infile3 = open('./use_file/Coupling_assay.pickle', 'rb')
coup = pickle.load(infile3)
rs3 = pd.concat([rs,rs1,rs2],axis=1)
co=[]
for index, row in rs3.iterrows():
    key = index.split('-')[0]  # extract the key from the index
    
    # check if the key exists in coup
    if key in coup:
        v = coup[key]
        co.append((index, row['Class'], row['GPCR_family'], row['gfam'], v[0], v[1], v[2], v[3]))
    else:
        # add a tuple with None values for the missing data
        co.append((index, row['Class'], row['GPCR_family'], row['gfam'], None, None, None, None))
coup1=pd.DataFrame(co,columns=['GPCR','Class','GPCR_family','gfam','Gs','Gio','Gq11','G1213']).set_index('GPCR')
rs3 = coup1
lut.update(lut1)
rs3 = rs3[ ['gfam'] + [ col for col in rs3.columns if col != 'gfam' ] ]
#coup={'G1213':'#800080','Only GtoPdb':'#cfccc6','Only TGF or GEMTA':'#FF8C00','Not coupled':'#808080'}
cou={'Only GtoPdb':'#cfccc6','Not coupled':'#808080'}
lut.update(cou) 

#####def clustermapp all
# Define the colormap for GPCR classes
def plot_clustermap(table,protein=str(),name=str()):    
    mask=0
    num_rows, num_cols = table.shape
    figsize = (num_cols //1.6, num_rows //1.6 )
    dendrogram_ratio = 0.15 * num_rows / num_cols
    cl = ListedColormap(['olive'])
    sns.set(style="ticks", rc={'axes.facecolor':'lightgray'})
    sns.set(font_scale=6.5)
    g3 = sns.clustermap(table, mask=table==mask, method='ward', row_colors=rs3, figsize=figsize, 
                        tree_kws=dict(linewidths=3, colors=(0, 0, 0)), dendrogram_ratio=(dendrogram_ratio, .1), 
                        colors_ratio= (0.01,0.01), col_cluster=False, cmap=cl, square=True, xticklabels=1, yticklabels=1,
                        linewidths=0.5, linecolor='silver')
    g3.ax_heatmap.set_xlabel('Position', fontsize=80, fontweight='bold')
    g3.ax_heatmap.set_ylabel('GPCR', fontsize=80, fontweight='bold')
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xmajorticklabels(), fontsize=80, fontweight='bold')
    g3.ax_heatmap.tick_params(axis='x', which='major', pad=2)
    g3.ax_heatmap.tick_params(axis='y', which='major', pad=2)
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_ymajorticklabels(), fontsize=80, fontweight='bold') 
    plt.setp(g3.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    g3.cax.set_visible(False)
    g3.gs.update(top=1.8, right=2)
    g3.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table.columns), linewidth=20, color='black')
    g3.ax_heatmap.hlines(y=len(table.index), xmin=-0.5, xmax=len(table.columns), linewidth=20, color='black')
    g3.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=20, color='black')
    g3.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=20, color='black')
    # gfamily and class

    handles = [Patch(facecolor=lut[name]) for name in lut]
    handles1 = [Patch(facecolor=lut2[name]) for name in lut2]
    handles2 = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
    handles3 = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
    hand=handles+handles2+handles3+handles1
    lt={**lut,  **gprot_ann, **ann_colors,**lut2}
    leg1 = g3.ax_heatmap.legend(hand, lt,frameon=False, ncol=6,loc='lower center',bbox_to_anchor=(1.1,-0.37), bbox_transform=plt.gcf().transFigure,fontsize=80, title='Gprotein family(UCM)/GPCR class/family/structure',title_fontsize=85, columnspacing=0.1, handletextpad=0.1, prop={'weight': 'bold'})
    # gpcr family
    if protein =='gprotein':
        for xtic in g3.ax_heatmap.get_xmajorticklabels():
            if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
                xtic.set_color(gp_ann[xtic.get_text()])    
    elif protein =='gpcr':
        for xtic in g3.ax_heatmap.get_xmajorticklabels():
            if xtic.get_text() in mydict.keys(): # Change color if exist else not
                xtic.set_color(mydict[xtic.get_text()])     
    elif protein =='pair':
        for xtic in g3.ax_heatmap.get_xmajorticklabels():
            if xtic.get_text().split('-')[0] in mydict.keys(): # Change color if exist else not
                xtic.set_color(mydict[xtic.get_text().split('-')[0]])            
    g3.savefig('./heatmaps/plots/Gs_vs_Gio_receptor_'+protein+'_cluster_'+name+'.svg',dpi=600,bbox_inches='tight',format='svg')


plot_clustermap(table1,'gprotein','all')
plot_clustermap(table12,'gpcr','all')
plot_clustermap(table1A,'gprotein','allA')
plot_clustermap(table12A,'gpcr','allA')
plot_clustermap(table12gp,'pair','02')
plot_clustermap(table1gAp,'pair','02_A')

#####def clustermapp 02
# Define the colormap for GPCR classes
def plot_clustermap1(table,protein=str(),name=str()):    
    mask=0
    num_rows, num_cols = table.shape
    figsize = (num_cols //1,num_rows //1.6 )
    dendrogram_ratio = 0.05 * num_rows / num_cols
    cl = ListedColormap(['olive'])
    sns.set(style="ticks", rc={'axes.facecolor':'lightgray'})
    sns.set(font_scale=6.5)
    g3 = sns.clustermap(table, mask=table==mask, method='ward', row_colors=rs3, figsize=figsize, 
                        tree_kws=dict(linewidths=3, colors=(0, 0, 0)), dendrogram_ratio=(dendrogram_ratio, .1), 
                        colors_ratio= (0.02,0.01), col_cluster=False, cmap=cl, square=True, xticklabels=1, yticklabels=1,
                        linewidths=0.5, linecolor='silver')
    g3.ax_heatmap.set_xlabel('Position', fontsize=80, fontweight='bold')
    g3.ax_heatmap.set_ylabel('GPCR', fontsize=80, fontweight='bold')
    g3.ax_heatmap.set_xticklabels(g3.ax_heatmap.get_xmajorticklabels(), fontsize=80, fontweight='bold')
    g3.ax_heatmap.tick_params(axis='x', which='major', pad=2)
    g3.ax_heatmap.tick_params(axis='y', which='major', pad=2)
    g3.ax_heatmap.set_yticklabels(g3.ax_heatmap.get_ymajorticklabels(), fontsize=80, fontweight='bold') 
    plt.setp(g3.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    g3.cax.set_visible(False)
    g3.gs.update(top=1.8, right=2)
    g3.ax_heatmap.hlines(y=0, xmin=-0.5, xmax=len(table.columns), linewidth=20, color='black')
    g3.ax_heatmap.hlines(y=len(table.index), xmin=-0.5, xmax=len(table.columns), linewidth=20, color='black')
    g3.ax_heatmap.vlines(x=0, ymin=0, ymax=len(table.index), linewidth=20, color='black')
    g3.ax_heatmap.vlines(x=len(table.columns), ymin=0, ymax=len(table.index), linewidth=20, color='black')
    # gfamily and class

    handles = [Patch(facecolor=lut[name]) for name in lut]
    handles1 = [Patch(facecolor=lut2[name]) for name in lut2]
    handles2 = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
    handles3 = [Patch(facecolor=gprot_ann[name]) for name in gprot_ann]
    hand=handles+handles2+handles3+handles1
    lt={**lut,  **gprot_ann, **ann_colors,**lut2}
    leg1 = g3.ax_heatmap.legend(hand, lt,frameon=False, ncol=2,loc='upper center',bbox_to_anchor=(2.9,1.6), bbox_transform=plt.gcf().transFigure,fontsize=80, title='Gprotein family(UCM)/GPCR class/family/structure',title_fontsize=85, columnspacing=0.1, handletextpad=0.1, prop={'weight': 'bold'})
    # gpcr family
    if protein =='gprotein':
        for xtic in g3.ax_heatmap.get_xmajorticklabels():
            if xtic.get_text() in gp_ann.keys(): # Change color if exist else not
                xtic.set_color(gp_ann[xtic.get_text()])    
    else:
        for xtic in g3.ax_heatmap.get_xmajorticklabels():
            if xtic.get_text() in mydict.keys(): # Change color if exist else not
                xtic.set_color(mydict[xtic.get_text()])     

    g3.savefig('./heatmaps/plots/Gs_vs_Gio_receptor_'+protein+'_cluster_'+name+'.svg',dpi=600,bbox_inches='tight',format='svg')

plot_clustermap1(table1g,'gprotein','02')
plot_clustermap1(table12g,'gpcr','02')
plot_clustermap1(table1gA,'gprotein','02A')
plot_clustermap1(table12gA,'gpcr','02A')
