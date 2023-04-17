
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from matplotlib_venn import venn3
import plotly.express as px
########## Create a sunburst diagram of the experimental set of structures annotating with class

fil=sys.argv[1]

df=pd.read_csv("./use_file/GPCR_structs_"+fil+".tsv", comment="#", sep="\t")
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein_Gene_name"].map(equiv)
#for venn
gpcr = df[["Receptor_Gene_name", "Gprotein_Gene_name","gfam"]].drop_duplicates()
gpcr1 = gpcr[gpcr.duplicated(subset=['Receptor_Gene_name'], keep=False)]
res = gpcr1.groupby('Receptor_Gene_name')['gfam'].apply(lambda x:'+'.join(x)).reset_index()
gpcr2 = pd.concat([gpcr, res], ignore_index=True)

#df2=pd.read_csv("targets_and_families.csv", comment="#", sep=",")


#gprotein based count
d=df.groupby(['gfam'])['PDB_ID'].count().reset_index()
d1=gpcr2.groupby(['gfam'])['Receptor_Gene_name'].count().reset_index()
d2=df.groupby(['Gprotein_Gene_name'])['Receptor_Gene_name'].nunique().reset_index()
d3=df.groupby(['Gprotein_Gene_name'])['PDB_ID'].nunique().reset_index()

#venn
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
color_purple='#800080'

gpcr2.set_index('gfam',inplace=True)
n1=gpcr2.loc['Gs', 'Receptor_Gene_name'].reset_index(drop=True)
n2=gpcr2.loc['Gio', 'Receptor_Gene_name'].reset_index(drop=True)
n3=gpcr2.loc['Gq11', 'Receptor_Gene_name'].reset_index(drop=True)
n4=gpcr2.loc['G1213', 'Receptor_Gene_name'].reset_index(drop=True)
d1=pd.concat([n1,n2,n3,n4],axis=1)

d1.set_axis(['Gs','Gio','Gq11','G1213'], axis="columns", inplace=True)
d1.to_csv("./statistics/sunburst/venn_gfam_"+fil+".tsv",index=None,header=True,sep='\t')
d3.to_csv("./statistics/sunburst/pdbid_gprotein_"+fil+".tsv",index=None,header=True,sep='\t')
d2.to_csv("./statistics/sunburst/gpcr_gprotein_"+fil+".tsv",index=None,header=True,sep='\t')
#class
d4=df.groupby(['Class'])['PDB_ID'].nunique().reset_index()
d5=df.groupby(['Class'])['Receptor_Gene_name'].nunique().reset_index()
d4.to_csv("./statistics/sunburst/pdbid_class_"+fil+".tsv",index=None,header=True,sep='\t')
d5.to_csv("./statistics/sunburst/gpcr_class_"+fil+".tsv",index=None,header=True,sep='\t')
d6=df.groupby(['Class','gfam'])['Receptor_Gene_name'].nunique().reset_index()
d6.to_csv("./statistics/sunburst/gpcr_class_gfam_"+fil+".tsv",index=None,header=True,sep='\t')
d7=df.groupby(['Class','gfam'])['PDB_ID'].nunique().reset_index()
d7.to_csv("./statistics/sunburst/pdbid_class_gfam_"+fil+".tsv",index=None,header=True,sep='\t')
d8=df.groupby(['Class','Gprotein_Gene_name'])['Receptor_Gene_name'].nunique().reset_index()
d8.to_csv("./statistics/sunburst/gprotein_gpcr_class_"+fil+".tsv",index=None,header=True,sep='\t')
d9=df.groupby(['Class','Gprotein_Gene_name'])['PDB_ID'].nunique().reset_index()
d9.to_csv("./statistics/sunburst/gprotein_pdbid_class_"+fil+".tsv",index=None,header=True,sep='\t')
d10=df.groupby(['Class','Gprotein_Gene_name','gfam'])['Receptor_Gene_name'].nunique().reset_index()
d11=df.groupby(['Class','Gprotein_Gene_name','gfam'])['PDB_ID'].nunique().reset_index()
d10.to_csv("./statistics/sunburst/gpcr_sunburst_"+fil+".tsv",index=None,header=True,sep='\t')
d11.to_csv("./statistics/sunburst/pdbid_sunburst_"+fil+".tsv",index=None,header=True,sep='\t')

import plotly.express as px
df=pd.read_csv("./statistics/sunburst/gpcr_sunburst_clean.tsv", comment="#", sep="\t")
lut = {'Gio':'#0000FF','Gq11':'#008000','Gs':'#FF0000','G1213':'#800080'}
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}

df["c1"] = df["gfam"].map(lut)
df["c2"] = df["Class"].map(lut1)
lut2=dict(zip(df.Gprotein_Gene_name, df.c1))
lut.update(lut1)
lut.update(lut2)

# create the nested donut chart using plotly
fig = px.sunburst(df, path=['gfam', 'Gprotein_Gene_name', 'Class'], values='Receptor_Gene_name',
                  color='gfam', color_discrete_map=lut, maxdepth=3,
                  branchvalues='total')
fig.update_traces(textinfo='label+value', texttemplate='%{label},%{value}' ,textfont_color='black',marker=dict(line=dict(width=0.5, color='white')))
fig.update_traces(marker_colors=[lut[cat] for cat in fig.data[-1].labels])
fig.update_layout(uniformtext=dict(minsize=16))
fig.show()
fig.write_image("./statistics/sunburst/gpcr_sunburst.svg", scale=6)


####
df=pd.read_csv("./statistics/sunburst/pdbid_sunburst_clean.tsv", comment="#", sep="\t")
lut = {'Gio':'#0000FF','Gq11':'#008000','Gs':'#FF0000','G1213':'#800080'}
lut1 = {'classA':'#44daeb','classB1':'#05fa98','classB2':'#0da813','classC':'#996f22','Frizzeled':'#ebd8b7'}

df["c1"] = df["gfam"].map(lut)
df["c2"] = df["Class"].map(lut1)
lut2=dict(zip(df.Gprotein_Gene_name, df.c1))
lut.update(lut1)
lut.update(lut2)

# create the nested donut chart using plotly
fig = px.sunburst(df, path=['gfam', 'Gprotein_Gene_name', 'Class'], values='PDB_ID',
                  color='gfam', color_discrete_map=lut, maxdepth=3,
                  branchvalues='total')
fig.update_traces(textinfo='label+value', texttemplate='%{label},%{value}' ,textfont_color='black',marker=dict(line=dict(width=0.5, color='white')))
fig.update_traces(marker_colors=[lut[cat] for cat in fig.data[-1].labels])
fig.update_layout(uniformtext=dict(minsize=16))
fig.show()
fig.write_image("./statistics/sunburst/pdbid_sunburst.svg", scale=6)



####gfam donut
df=pd.read_csv("./use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein_Gene_name"].map(equiv)
#for venn
gpcr = df[["Receptor_Gene_name", "Gprotein_Gene_name","gfam"]].drop_duplicates()
gpcr1 = gpcr[gpcr.duplicated(subset=['Receptor_Gene_name'], keep=False)]
d=gpcr.groupby(['gfam'])['Receptor_Gene_name'].count().reset_index()

color_red = '#FF0000'
color_blue = '#0000FF'
color_green='#008000'
color_purple='#800080'
palette_color = [color_red,color_blue,color_green,color_purple]
sizes = list(d.Receptor_Gene_name)
labels = list(d.gfam)
colors =[color_purple,color_blue,color_green,color_red]
plt.figure(figsize=(30,30))

p, tx, autotexts =plt.pie(sizes, colors=colors, labels=labels, pctdistance=0.85,
        autopct="",textprops={'fontsize': 90},startangle=45)
for w in p:
    w.set_linewidth(2)
    w.set_edgecolor('black')
centre_circle = plt.Circle((0, 0), 0.6, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes[i]))

plt.savefig('./statistics/sunburst/gfam_gpcr_pie.svg',bbox_inches='tight', dpi=600, format='svg')



