#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd


df=pd.read_csv("./use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
df1=df.groupby(['Receptor_Gene_name','Gprotein_Gene_name'])['Resolution'].min().reset_index()
df2= df1.merge(df, on=['Receptor_Gene_name','Gprotein_Gene_name','Resolution'],how="left")
df3=df2.sort_values('Gprotein_coverage',ascending=False).drop_duplicates(subset=['Receptor_Gene_name','Gprotein_Gene_name'], keep='first')
c=list(df)
df3=df3[list(df)]
df3.to_csv("./use_file/GPCR_structs_clean_consensus.tsv",index=None,header=True,sep='\t')

