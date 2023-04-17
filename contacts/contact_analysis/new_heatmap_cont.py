from matplotlib import pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sys
import pandas as pd
import numpy as np
#### Heatmap of most conserved contacts in alphafold analysis for all Gprotein family and comparison of contacts with the ones in experimental structures
fil='all'
df=pd.read_csv("gio_contact_"+fil+".tsv", comment="#", sep="\t")
df = df[(df['unique_gpcr_contact']/len(set(df['GPCR'])) > 0.2)]
df1=pd.read_csv("gs_contact_"+fil+".tsv", comment="#", sep="\t")
df1 = df1[(df1['unique_gpcr_contact']/len(set(df1['GPCR'])) > 0.2)]
df2=pd.read_csv("g1213_contact_"+fil+".tsv", comment="#", sep="\t")
df2 = df2[(df2['unique_gpcr_contact']/len(set(df2['GPCR'])) > 0.2)]

df3=pd.read_csv("gq11_contact_"+fil+".tsv", comment="#", sep="\t")
df3 = df3[(df3['unique_gpcr_contact']/len(set(df3['GPCR'])) > 0.2)]

d_all=pd.concat([df,df1,df2,df3]).drop_duplicates(keep=False)
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
t=d_all[["BW","hex"]].drop_duplicates()
mydict = dict(zip(t.BW, t.hex))


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
g1213=d_all[d_all.gfam =='G1213']
gq11=d_all[d_all.gfam =='Gq11']

#gio
pos=d_all[['BW','Gprot_pos','Structural','Gprot_struct']].drop_duplicates()
gio = pos.merge(gio,on=['BW','Gprot_pos','Structural','Gprot_struct'],how="left")
sorter1=['HN','hns1','S1','s1h1','H','h1ha','HA','hahb','HB','hbhc','HC','hchd','HD','hdhe','HE','hehf','HF','hfs2','S2','s2s3','S3','s3h2','H2','h2s4','S4','s4h3','H3','h3s5','S5','s5hg','HG','hgh4','H4','h4s6','S6','s6h5','H5']
gio.Gprot_struct = gio.Gprot_struct.astype("category")
gio.Gprot_struct.cat.set_categories(sorter1, inplace=True)
gio[['G','num']] = gio['Gprot_pos'].str.split('.', expand=True)
gio['num']=gio['num'].astype(float)
ds1=gio[['Gprot_struct','Gprot_pos','num']].drop_duplicates()
ds1.sort_values(by=['Gprot_struct',"num"],ascending=True,inplace=True)
ls2=ds1["Gprot_pos"].drop_duplicates().tolist()
sorter=['Nterm','TM1','ICL1','TM2','ECL1','TM3','ICL2','TM4','ECL2','TM5','ICL3','TM6','ECL3','TM7','H8','Cterm']
gio.Structural = gio.Structural.astype("category")
gio.Structural.cat.set_categories(sorter, inplace=True)
gio.sort_values(["Structural","BW"],inplace=True)
ls1=gio["BW"].drop_duplicates().tolist()
gio['Gprot_pos'] = pd.Categorical(gio['Gprot_pos'],
                                   categories=ls2,
                                   ordered=True)

def triangulation_for_triheatmap(M, N):
    xv, yv = np.meshgrid(np.arange(-0.5, M), np.arange(-0.5, N))  # vertices of the little squares
    xc, yc = np.meshgrid(np.arange(0, M), np.arange(0, N))  # centers of the little squares
    x = np.concatenate([xv.ravel(), xc.ravel()])
    y = np.concatenate([yv.ravel(), yc.ravel()])
    cstart = (M + 1) * (N + 1)  # indices of the centers

    trianglesN = [(i + j * (M + 1), i + 1 + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesE = [(i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesS = [(i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    trianglesW = [(i + (j + 1) * (M + 1), i + j * (M + 1), cstart + i + j * M)
                  for j in range(N) for i in range(M)]
    return [Triangulation(x, y, triangles) for triangles in [trianglesN, trianglesE, trianglesS, trianglesW]]

tab = pd.pivot_table(d_all, values='fraction_pair', index=['Gprot_pos','BW'],columns=['gfam']).reset_index()
df_piv = tab.pivot_table(index='Gprot_pos', columns='BW', dropna=False)
df_piv= df_piv.reindex(ls2)
l1=[]
for i in ls1:
    for j in ['G1213', 'Gio', 'Gq11', 'Gs']:
        l1.append((j,i))

df_piv = df_piv.reindex(columns=l1)
M = len(df_piv.columns) // 4
N = len(df_piv)
values = [df_piv[dir] for dir in
          ['G1213', 'Gio', 'Gq11', 'Gs']]  # these are the 4 column names in df
triangul = triangulation_for_triheatmap(M, N)
cmap0 = plt.cm.get_cmap("Purples").copy()
cmap1 = plt.cm.get_cmap("Blues").copy()
cmap2 = plt.cm.get_cmap("Greens").copy()
cmap3 = plt.cm.get_cmap("Reds").copy()
cmaps = [cmap0,cmap1,cmap2,cmap3] 

# Choose the color
cmap0.set_bad(color='snow')
cmap1.set_bad(color='snow')
cmap2.set_bad(color='snow')
cmap3.set_bad(color='snow')

fig, ax = plt.subplots(figsize=(55, 55))
imgs = [ax.tripcolor(t, np.ravel(val), cmap=cmap, vmin=0, vmax=1, ec='white',lw=0.00000001)
        for t, val, cmap in zip(triangul, values, cmaps)]
ax.tick_params(length=1)
ax.set_xticks(range(M))
ax.set_xticklabels(df_piv['G1213'].columns, fontsize=40,rotation=90, fontweight='bold')
ax.set_yticks(range(N))
ax.set_yticklabels(df_piv.index, fontsize=40, fontweight='bold')
ax.margins(x=0, y=0)
ax.set_aspect('equal', 'box')  # square cells

[plt.axhline(y=i, linestyle='-',linewidth=1,color='black') for i in np.arange(-0.5,N,1)]
[plt.axvline(x=i, linestyle='-',linewidth=1,color='black') for i in np.arange(-0.5,M,1)]

for ytic in ax.get_yticklabels():
    if ytic.get_text() in gp_ann.keys(): # Change color if exist else not
        ytic.set_color(gp_ann[ytic.get_text()])
for xtic in ax.get_xticklabels():
    if xtic.get_text() in mydict.keys(): # Change color if exist else not
        xtic.set_color(mydict[xtic.get_text()])
ann_colors.update(gprot_ann) 
handles = [Patch(facecolor=ann_colors[name]) for name in ann_colors]
leg2=ax.legend(handles, ann_colors,bbox_to_anchor=(0.52,0.68), bbox_transform=plt.gcf().transFigure,frameon=False, loc='best',ncol=1,fontsize=50,title='GPCR_pos/Gprot_pos',title_fontsize=50,columnspacing=0.04,labelspacing=0.01,handletextpad =0.01)


cbar=plt.colorbar(imgs[0], pad=0.002)
cbar.set_ticks(np.arange(0, 1.1, 0.1))
cbar.ax.set_position([0.055,0.45, 0.0035, 0.2])
cbar.ax.tick_params(labelsize=35)
cbar1=plt.colorbar(imgs[1], pad=0.002)
cbar1.set_ticks(np.arange(0, 1.1, 0.1))
cbar1.ax.set_position([0.065,0.45, 0.0035, 0.2])
cbar1.ax.tick_params(labelsize=35)
cbar2=plt.colorbar(imgs[2], pad=0.002)
cbar2.set_ticks(np.arange(0, 1.1, 0.1))
cbar2.ax.set_position([0.045,0.45, 0.0035, 0.2])
cbar2.ax.tick_params(labelsize=35)
cbar3=plt.colorbar(imgs[3], pad=0.002)
cbar3.set_ticks(np.arange(0, 1.1, 0.1))
cbar3.ax.set_position([0.035,0.45, 0.0035, 0.2])
cbar3.ax.tick_params(labelsize=35)
cbar.ax.set_title('G1213',fontsize=40,pad=20,rotation=90)
cbar1.ax.set_title('Gio',fontsize=40,pad=20,rotation=90)
cbar2.ax.set_title('Gq11',fontsize=40,pad=20,rotation=90)
cbar3.ax.set_title('Gs',fontsize=40,pad=20,rotation=90)
cbar.set_ticks([])
cbar2.set_ticks([])
cbar3.set_ticks([])

fp=pd.read_csv("fraction_pair_all.tsv", comment="#", sep="\t")
fp = fp.drop(fp[fp.fraction_pair < 0.2].index)
gpc={k: v for v, k in enumerate(ls1)}
gpr={k: v for v, k in enumerate(ls2)}
fp["gpr"] = fp["Gprot_pos"].map(gpr)
fp["gpc"] = fp["BW"].map(gpc)
fp=fp.dropna()



###G1213
g12=tab[["Gprot_pos","BW","G1213"]].drop_duplicates().dropna()
g1213 = fp[fp['gfam'] =='G1213']
g121=g12.merge(g1213,left_on=['BW','Gprot_pos'],right_on=['BW','Gprot_pos'],how="left").dropna()

for dir in [(-1, 0)]:
    for i,j in g121.iterrows():
        ax.text(j['gpc'] + 0.3 * dir[1], j['gpr'] + 0.45 * dir[0],"*", color='k', ha='center', va='center',fontsize=22, weight='bold')
##Gio
###G1213
gi=tab[["Gprot_pos","BW","Gio"]].drop_duplicates().dropna()
gio = fp[fp['gfam'] =='Gio']
gio1=gi.merge(gio,left_on=['BW','Gprot_pos'],right_on=['BW','Gprot_pos'],how="left").dropna()
for dir in [(0, 1)]:
    for i,j in gio1.iterrows():
        ax.text(j['gpc'] + 0.3 * dir[1], j['gpr'] -0.5 * dir[0],"*", color='k', ha='center', va='center',fontsize=22, weight='bold')

##Gq11
gq=tab[["Gprot_pos","BW","Gq11"]].drop_duplicates().dropna()
gq1 = fp[fp['gfam'] =='Gq11']
gq11=gq.merge(gq1,left_on=['BW','Gprot_pos'],right_on=['BW','Gprot_pos'],how="left").dropna()
for dir in [(1, 0)]:
    for i,j in gq11.iterrows():
        ax.text(j['gpc'] + 0.3 * dir[1], j['gpr'] + 0.15 * dir[0],"*", color='k', ha='center', va='center',fontsize=22, weight='bold')
##Gs
gs=tab[["Gprot_pos","BW","Gs"]].drop_duplicates().dropna()
gs1 = fp[fp['gfam'] =='Gs']
gs11=gs.merge(gs1,left_on=['BW','Gprot_pos'],right_on=['BW','Gprot_pos'],how="left").dropna()
for dir in [(0, -1)]:
    for i,j in gs11.iterrows():
        ax.text(j['gpc'] + 0.33 * dir[1], j['gpr'] -0.5 * dir[0],"*", color='k', ha='center', va='center',fontsize=22, weight='bold')









fig.savefig('heatmap_all.svg',dpi=600,bbox_inches='tight',format='svg')
