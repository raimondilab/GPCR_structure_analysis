from scipy.cluster import hierarchy
from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
from matplotlib.patches import Patch
import seaborn as sns; sns.set(color_codes=True)
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import matplotlib.patches as mpatches

from mpl_toolkits.axes_grid1 import make_axes_locatable
from statannot import add_stat_annotation
import sys 
import numpy as np
import pickle
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
ab=ab.fillna('#da70d6')

ab1={'Antibody':'#da70d6','No antibody':'#f2f2f2'}
lut.update(ab1) 
ab=ab.drop(['Antibody'], axis=1)
rs3 = pd.concat([rs3,ab],axis=1)




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
num_rows, num_cols = table1.shape
figsize = (num_cols //4, num_rows //4)
mask=0
g3=sns.clustermap(table1,method='ward',mask=table1==mask,col_colors=rs3,row_colors=rs3,tree_kws=dict(linewidths=2, colors=(0, 0, 0)),dendrogram_ratio=(.08, 0.08),cmap='Reds',square=True,xticklabels=1,yticklabels=1,linewidths=0.5,linecolor='silver',colors_ratio=(0.007,0.007),figsize=figsize)

linkage_matrix = g3.dendrogram_row.linkage
clusters = hierarchy.fcluster(linkage_matrix, 3, criterion='maxclust')
clust={}
for idx, x in enumerate(clusters):
    for k,v in enumerate(table1.index.tolist()):
        if idx == k:
            clust[v]=x
rs4=rs3.copy()        
rs4["clust"] = rs4.index.map(clust)
rs4=rs4.drop(columns=['Class', 'antibody'])
inv_map = {v: k for k, v in lut.items()}

rs4.replace(inv_map, inplace=True)
cl1 = rs4[rs4['clust'] == 1]
cl2 = rs4[rs4['clust'] == 2]
cl3 = rs4[rs4['clust'] == 3]
ccl1=cl1.apply(pd.Series.value_counts).fillna(0)
ccl2=cl2.apply(pd.Series.value_counts).fillna(0)
ccl3=cl3.apply(pd.Series.value_counts).fillna(0)

# get clusters1
Gs1=(ccl1["Gs"].loc[ 'Only GtoPdb']+ccl1["Gs"].loc[ 'Gs'])/ccl1['Gs'].sum()
Gio1=(ccl1["Gio"].loc[ 'Only GtoPdb']+ccl1["Gio"].loc[ 'Gio'])/ccl1['Gs'].sum()
Gq111=(ccl1["Gq11"].loc[ 'Only GtoPdb']+ccl1["Gq11"].loc[ 'Gq11'])/ccl1['Gs'].sum()
G12131=(ccl1["G1213"].loc[ 'Only GtoPdb']+ccl1["G1213"].loc[ 'G1213'])/ccl1['Gs'].sum()
####
Gs2=(ccl2["Gs"].loc[ 'Only GtoPdb']+ccl2["Gs"].loc[ 'Gs'])/ccl2['Gs'].sum()
Gio2=(ccl2["Gio"].loc[ 'Only GtoPdb']+ccl2["Gio"].loc[ 'Gio'])/ccl2['Gs'].sum()
Gq112=(ccl2["Gq11"].loc[ 'Only GtoPdb']+ccl2["Gq11"].loc[ 'Gq11'])/ccl2['Gs'].sum()
G12132=(ccl2["G1213"].loc[ 'Only GtoPdb'])/ccl2['Gs'].sum()
####
Gs3=(ccl3["Gs"].loc[ 'Only GtoPdb']+ccl3["Gs"].loc[ 'Gs'])/ccl3['Gs'].sum()
Gio3=(ccl3["Gio"].loc[ 'Only GtoPdb']+ccl3["Gio"].loc[ 'Gio'])/ccl3['Gs'].sum()
Gq113=(ccl3["Gq11"].loc[ 'Only GtoPdb']+ccl3["Gq11"].loc[ 'Gq11'])/ccl3['Gs'].sum()
G12133=(ccl3["G1213"].loc[ 'Only GtoPdb']+ccl3["G1213"].loc[ 'G1213'])/ccl3['Gs'].sum()



ratio = pd.DataFrame({
    'gfam': ['Gs', 'Gio', 'Gq11', 'G1213','Gs', 'Gio', 'Gq11', 'G1213','Gs', 'Gio', 'Gq11', 'G1213'],
    'val': [Gs1, Gio1, Gq111, G12131,Gs2, Gio2, Gq112, G12132,Gs3, Gio3, Gq113, G12133],
    'cluster':['cluster1','cluster1','cluster1','cluster1','cluster2','cluster2','cluster2','cluster2','cluster3','cluster3','cluster3','cluster3']})
#####plot
ratio.to_csv('./plots/ratio_cluster_rmsd.tsv',index=None,sep='\t')


color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
color_purple='#800080'
colors = [color_red,color_blue,color_green,color_purple,color_red,color_blue,color_green,color_purple,color_red,color_blue,color_green,color_purple]

#cluster1
sns.set(style='ticks')
plt.figure(figsize=(15,15))

g = sns.barplot(
    data=ratio,
    x="gfam", y="val", alpha=1,palette=colors,hue='cluster', edgecolor = 'black')

plt.xlabel("Gprotein Family",fontsize=60)
plt.ylabel("Fraction of coupling",fontsize=60)
plt.xticks(fontsize=60,rotation=90)
plt.yticks(np.arange(0,1.1,0.1),fontsize=60)
plt.title('Clusters',fontsize=60)
hatches = ['\\','\\','\\','\\','-','-','-','-','x','x','x','x']
# Loop over the bars
for i,thisbar in enumerate(g.patches):
    # Set a different hatch for each bar
    thisbar.set_color(colors[i])
for i,thisbar in enumerate(g.patches):
    thisbar.set_edgecolor("black")
    thisbar.set_linewidth(3)
# Loop over the bars
for i,s in enumerate(g.patches):
    # Set a different hatch for each bar
    s.set_hatch(hatches[i])
clu1 = mpatches.Patch(edgecolor='black', facecolor='white',
                     hatch='\\', label='Cluster 1')
clu2 = mpatches.Patch(edgecolor='black', facecolor='white', label='Cluster 2',hatch='+')
clu3 = mpatches.Patch(edgecolor='black', facecolor='white', label='Cluster 3',hatch='x')

plt.legend(handles=[clu1,clu2,clu3],bbox_to_anchor=(0.66,0.98),fontsize=30,title='Gprotein Family',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0.1 ,framealpha=1)
plt.savefig('./plots/clusters_all_ratio.pdf',bbox_inches='tight', dpi=600, format='pdf')
