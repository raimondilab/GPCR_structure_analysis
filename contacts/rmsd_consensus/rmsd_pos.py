#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
import numpy as np
df=pd.read_csv("GPCR_structs_clean_consensus.tsv", comment="#", sep="\t")
gprot=pd.read_csv("./use_file/human_sequence_table_FR_new.txt", comment="#", sep="\t",index_col=0)
gprot1=gprot.unstack().reset_index()
gprot1.columns=['Gprotein_id','CGN','Pos2']
gprot1['Pos2'] = gprot1['Pos2'].str.replace(r'\D', '')
gprot1['Pos2'] .replace('', np.nan, inplace=True)
gprot1=gprot1.dropna()
ds=df[['Gprotein_Uniprot_AC', 'Gprotein_Gene_name','Gprotein_Uniprot_ID']].drop_duplicates()
gprot2 = gprot1.merge(ds,left_on=['Gprotein_id'],right_on=['Gprotein_Uniprot_ID'],how="left")
gprot2=gprot2.dropna()
gprot2['Pos2']= gprot2['Pos2'].astype(int)
df2=pd.read_csv("./use_file/output_gpcrdb.tsv", comment="#", sep="\t")
df3=pd.read_csv("uniprot_pdb.tsv", sep="\t")
df3['UNIPROT_RES_NUM']= df3['UNIPROT_RES_NUM'].astype(int)
df3['CHAIN'][(df3['#PDBID']=='7e9h')&(df3['CHAIN']=='S')] ='R'
df3['CHAIN'][(df3['#PDBID']=='7e9g')&(df3['CHAIN']=='S')] ='R'
df3['CHAIN'][(df3['#PDBID']=='7mts')&(df3['CHAIN']=='B')] ='A'

df4 = df3.merge(gprot2,left_on=['UNIPROT_AC','UNIPROT_RES_NUM'],right_on=['Gprotein_Uniprot_AC','Pos2'],how="left")
df4=df4.dropna()
df4 = df4.loc[~((df4['#PDBID'] == '7mby') & (df4['Gprotein_Gene_name'] == 'GNAI1')),:]
d = df.merge(df3,right_on=['#PDBID','CHAIN'],left_on=['PDB_ID','Receptor_Chain'],how="left")
d=d.dropna()
df5 = d.merge(df2,left_on=['UNIPROT_AC','UNIPROT_RES_NUM'],right_on=['Uniprot','Pos1'],how="left")
df5=df5.dropna()
df5=df5.drop_duplicates()
df6=df5.loc[df5['BW'].isin(['1.49','1.5','1.51','3.49','3.5','3.51','4.49','4.5','4.51','5.41','6.49','6.5','6.51'])]
#gprot all
d6=df4.groupby(['CGN'])["#PDBID"].count().reset_index()
d6=d6.loc[d6['#PDBID'].isin([df['PDB_ID'].unique().size])]
d6.to_csv('gprot_pos_set_all.tsv',index=None,header=True,sep='\t')

df7=df4.loc[df4['CGN'].isin(d6['CGN'].tolist())]
#gpcr all
d5=df5.groupby(['BW'])['PDB_ID'].count().reset_index()
d5=d5.loc[d5['PDB_ID'].isin([df['PDB_ID'].unique().size])]
d5.to_csv('gpcr_pos_set_all.tsv',index=None,header=True,sep='\t')
df8=df5.loc[df5['BW'].isin(d5['BW'].tolist())]
df6['PDB_RES_NUM']= df6['PDB_RES_NUM'].astype(str)
df7['PDB_RES_NUM']= df7['PDB_RES_NUM'].astype(str)
df8['PDB_RES_NUM']= df8['PDB_RES_NUM'].astype(str)
ds6=df6.groupby(['PDB_ID','CHAIN'])["PDB_RES_NUM"].apply(lambda item:','.join(item)).reset_index()
ds7=df7.groupby(['#PDBID','CHAIN'])["PDB_RES_NUM"].apply(lambda item:','.join(item)).reset_index()
ds8=df8.groupby(['PDB_ID','CHAIN'])["PDB_RES_NUM"].apply(lambda item:','.join(item)).reset_index()
ds6.to_csv('gpcr_pos_cons_12.tsv',index=None,header=True,sep='\t')
ds7.to_csv('gprot_pos_cons.tsv',index=None,header=True,sep='\t')
ds8.to_csv('gpcr_pos_cons.tsv',index=None,header=True,sep='\t')
