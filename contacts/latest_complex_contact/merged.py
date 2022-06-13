#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
from functools import reduce
df1=pd.read_csv('./network/files/gs_nodes_all.tsv', comment="#", sep="\t")
df2=pd.read_csv('./network/files/gio_nodes_all.tsv', comment="#", sep="\t")
df3=pd.read_csv('./network/files/gs_nodes_all_classA.tsv', comment="#", sep="\t")
df4=pd.read_csv('./network/files/gio_nodes_all_classA.tsv', comment="#", sep="\t")

df5=pd.read_csv('./network/files/gq11_nodes_all.tsv', comment="#", sep="\t")

#df6=pd.read_csv('./network/files/g1213_nodes_all.tsv', comment="#", sep="\t")
df7=pd.read_csv('./network/files/gq11_nodes_all_classA.tsv', comment="#", sep="\t")
#df8=pd.read_csv(./network/files/'g1213_nodes_all_classA.tsv', comment="#", sep="\t")




data_frames = [df1, df2,df5]

df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['struct'],
                                            how='outer'), data_frames)
df_merged=df_merged.fillna(0)
df_merged.to_csv('./network/files/Nodes_struct_all.tsv',index=None,header=True,sep='\t')

d1=pd.read_csv('./network/files/gs_links_all.tsv', comment="#", sep="\t")
d2=pd.read_csv('./network/files/gio_links_all.tsv', comment="#", sep="\t")
d3=pd.read_csv('./network/files/gq11_links_all.tsv', comment="#", sep="\t")
#d4=pd.read_csv('./network/files/g1213_links_all.tsv', comment="#", sep="\t")

data_frames1 = [d1, d2,d3]

df_merged1 = reduce(lambda  left,right: pd.merge(left,right,on=['SS1','SS2'],
                                            how='outer'), data_frames1)
df_merged1=df_merged1.fillna(0)
df_merged1.to_csv('./network/files/Links_struct_all.tsv',index=None,header=True,sep='\t')


#classA
data_frames2 = [df3, df4,df7]
df_merged2 = reduce(lambda  left,right: pd.merge(left,right,on=['struct'],how='outer'), data_frames2)
df_merged2=df_merged2.fillna(0)
df_merged2.to_csv('./network/files/Nodes_struct_all_classA.tsv',index=None,header=True,sep='\t')
d5=pd.read_csv('./network/files/gs_links_all_classA.tsv', comment="#", sep="\t")
d6=pd.read_csv('./network/files/gio_links_all_classA.tsv', comment="#", sep="\t")

d7=pd.read_csv('./network/files/gq11_links_all_classA.tsv', comment="#", sep="\t")
#d8=pd.read_csv('./network/files/g1213_links_all_classA.tsv', comment="#", sep="\t")

data_frames3 = [d5,d6,d7]
df_merged3 = reduce(lambda  left,right: pd.merge(left,right,on=['SS1','SS2'],how='outer'), data_frames3)
df_merged3=df_merged3.fillna(0)
df_merged3.to_csv('./network/files/Links_struct_all_classA.tsv',index=None,header=True,sep='\t')

