#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import seaborn as sns
from matplotlib_venn import venn3


fil=sys.argv[1]

df=pd.read_csv("./use_file/GPCR_structs_"+fil+".tsv", comment="#", sep="\t")
equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein_Gene_name"].map(equiv)
df1=pd.read_csv("./use_file/classification.txt", comment="#", sep="\t")
#for venn
gpcr = df[["Receptor_Gene_name", "gfam"]].drop_duplicates()
gpcr1 = gpcr[gpcr.duplicated(subset=['Receptor_Gene_name'], keep=False)]



res = gpcr1.groupby('Receptor_Gene_name')['gfam'].apply(lambda x:'+'.join(x)).reset_index()
gpcr2 = pd.concat([gpcr, res], ignore_index=True)

#df2=pd.read_csv("targets_and_families.csv", comment="#", sep=",")


#gprotein based count
d=df.groupby(['gfam'])['PDB_ID'].count().reset_index()
d1=gpcr2.groupby(['gfam'])['Receptor_Gene_name'].nunique().reset_index()
d2=df.groupby(['Gprotein_Gene_name'])['Receptor_Gene_name'].nunique().reset_index()
d3=df.groupby(['Gprotein_Gene_name'])['PDB_ID'].nunique().reset_index()

#venn
color_red = '#c24a4d'
color_blue = '#4c72b0'
color_green='#008000'
d1.set_index('gfam',inplace=True)
n1=d1.loc['Gs', 'Receptor_Gene_name'] 
n2=d1.loc['Gio', 'Receptor_Gene_name']
n3=d1.loc['Gq11', 'Receptor_Gene_name']
n4=d1.loc['Gio+Gs', 'Receptor_Gene_name']+d1.loc['Gs+Gio', 'Receptor_Gene_name']
n5=d1.loc['Gio+Gq11', 'Receptor_Gene_name']
n6=d1.loc['Gio+Gs+Gq11', 'Receptor_Gene_name']


ax=plt.figure(figsize=(10,10))

ax=venn3(subsets = (n2,n1,n4,n3,n5,n6,1), set_labels = ('Gio', 'Gs','Gq11'),set_colors=(color_blue,color_red,color_green ))
for t in ax.set_labels: t.set_fontsize(50)
for t in ax.subset_labels: t.set_fontsize(30)
plt.savefig('./statistics/plots/gpcr_gfam_venn_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

#pie


sizes = list(d.PDB_ID)
labels = list(d.gfam)
colors = [color_blue,color_green,color_red]
explode = (0.05, 0.05,0.05)
ax=plt.figure(figsize=(10,10))

p, tx, autotexts =plt.pie(sizes, colors=colors, labels=labels, pctdistance=0.85,
        explode=explode, autopct="",textprops={'fontsize': 45})

centre_circle = plt.Circle((0, 0), 0.70, fc='white')
fig = plt.gcf()
fig.gca().add_artist(centre_circle)

for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes[i]))
plt.savefig('./statistics/plots/pdbid_gfam_pie_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


#pie2
palette_color = sns.color_palette('tab10')
sizes = list(d3.PDB_ID)
labels = list(d3.Gprotein_Gene_name)
colors = [color_blue,color_green,color_red]
plt.figure(figsize=(30,30))

p, tx, autotexts =plt.pie(sizes, colors=palette_color, labels=labels, pctdistance=2,
        autopct="",textprops={'fontsize': 37})
for w in p:
    w.set_linewidth(2)
    w.set_edgecolor('black')
centre_circle = plt.Circle((0, 0), 0.1, fc='white')
fig = plt.gcf()
plt.savefig('./statistics/plots/pdbid_gprotein_pie_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


#pie3
palette_color = sns.color_palette('tab10')
sizes = list(d2.Receptor_Gene_name)
labels = list(d2.Gprotein_Gene_name)
colors = [color_blue,color_red]
plt.figure(figsize=(30,30))

p, tx, autotexts =plt.pie(sizes, colors=palette_color, labels=labels, pctdistance=0.7,
        autopct="",textprops={'fontsize': 37})
for w in p:
    w.set_linewidth(2)
    w.set_edgecolor('black')
centre_circle = plt.Circle((0, 0), 0.1, fc='white')
fig = plt.gcf()
plt.savefig('./statistics/plots/gpcr_gprotein_pie_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




#class
df3= df.merge(df1, left_on=['Receptor_Gene_name'],right_on=['GPCR'],how="left")
d4=df3.groupby(['class'])['PDB_ID'].nunique().reset_index()
d5=df3.groupby(['class'])['GPCR'].nunique().reset_index()




#pieclass
palette_color = sns.color_palette('tab10')
sizes = list(d4.PDB_ID)
labels = list(d4['class'])
plt.figure(figsize=(10,10))

p, tx, autotexts =plt.pie(sizes, colors=palette_color, labels=labels, pctdistance=0.7,
        autopct="",textprops={'fontsize': 30})
for w in p:
    w.set_linewidth(2)
    w.set_edgecolor('black')
centre_circle = plt.Circle((0, 0), 0.1, fc='white')
fig = plt.gcf()
plt.savefig('./statistics/plots/pdbid_class_pie_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

#pieclass
palette_color = sns.color_palette('tab10')
sizes = list(d5.GPCR)
labels = list(d5['class'])
plt.figure(figsize=(10,10))

p, tx, autotexts =plt.pie(sizes, colors=palette_color, labels=labels, pctdistance=0.7,
        autopct="",textprops={'fontsize': 30})
for w in p:
    w.set_linewidth(2)
    w.set_edgecolor('black')
centre_circle = plt.Circle((0, 0), 0.1, fc='white')
fig = plt.gcf()
plt.savefig('./statistics/plots/gpcr_class_pie_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')




d6=df3.groupby(['class','gfam'])['GPCR'].nunique().reset_index()
d6=d6.sort_values(['GPCR'],ascending=False).reset_index(drop=True)
palette ={"classA": "C0", "classB1": "C1", 'classB2': "C2", "classC": "C3", "Frizzeled": "C4"}

# Draw a nested barplot by species and sex
sns.set(style='ticks')
plt.figure(figsize=(10,10))

g = sns.barplot(
    data=d6,
    x="gfam", y="GPCR", hue='class'
   , palette=palette, alpha=1)
plt.xlabel("Gprotein Family",fontsize=50)
plt.ylabel("Number of GPCR",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(1.42,1),fontsize=30,title='GPCR Class',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)

plt.savefig('./statistics/plots/gpcr_class_pie_bar_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


d7=df3.groupby(['class','gfam'])['PDB_ID'].nunique().reset_index()

d7=d7.sort_values(['PDB_ID'],ascending=False).reset_index(drop=True)
palette ={"classA": "C0", "classB1": "C1", 'classB2': "C2", "classC": "C3", "Frizzeled": "C4"}

sns.set(style='ticks')
plt.figure(figsize=(10,10))

g = sns.barplot(
    data=d7,
    x="gfam", y="PDB_ID", hue="class"
   , palette=palette, alpha=1)
plt.xlabel("Gprotein Family",fontsize=50)
plt.ylabel("Number of Structures",fontsize=50)
plt.xticks(fontsize=50)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(1.42,1),fontsize=30,title='GPCR Class',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./statistics/plots/pdbid_class_pie_bar_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')

d8=df3.groupby(['class','Gprotein_Gene_name'])['GPCR'].nunique().reset_index()
d8=d8.sort_values(['GPCR'],ascending=False).reset_index(drop=True)

palette ={"classA": "C0", "classB1": "C1", 'classB2': "C2", "classC": "C3", "Frizzeled": "C4"}

# Draw a nested barplot by species and sex
sns.set(style='ticks')
plt.figure(figsize=(10,10))

g = sns.barplot(
    data=d8,
    x="Gprotein_Gene_name", y="GPCR", hue="class"
   , palette=palette, alpha=1)
plt.xlabel("Gprotein",fontsize=50)
plt.ylabel("Number of GPCR",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(1.42,1),fontsize=30,title='GPCR Class',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./statistics/plots/gprotein_gpcr_class_pie_bar_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')


d9=df3.groupby(['class','Gprotein_Gene_name'])['PDB_ID'].nunique().reset_index()
d9=d9.sort_values(['PDB_ID'],ascending=False).reset_index(drop=True)

palette ={"classA": "C0", "classB1": "C1", 'classB2': "C2", "classC": "C3", "Frizzeled": "C4"}

# Draw a nested barplot by species and sex
sns.set(style='ticks')
plt.figure(figsize=(10,10))

g = sns.barplot(
    data=d9,
    x="Gprotein_Gene_name", y="PDB_ID", hue="class"
   , palette=palette, alpha=1)
plt.xlabel("Gprotein",fontsize=50)
plt.ylabel("Number of structures",fontsize=50)
plt.xticks(fontsize=50,rotation=90)
plt.yticks(fontsize=50)
plt.legend(bbox_to_anchor=(1.42,1),fontsize=30,title='GPCR Class',title_fontsize=35,ncol=1,borderpad=0.2,edgecolor='black', borderaxespad=0 ,framealpha=1)
plt.savefig('./statistics/plots/gprotein_pdbid_class_pie_bar_'+fil+'.svg',bbox_inches='tight', dpi=600, format='svg')





