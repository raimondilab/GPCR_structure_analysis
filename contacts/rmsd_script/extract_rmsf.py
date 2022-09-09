#!/data/SW/anaconda3/envs/myenv/bin/python
import sys, os, operator,gzip,glob, math
import pandas as pd
path=sys.argv[1]
all_files = glob.glob("/data/Users/marin/contact_analysis/GPCR_structure_analysis/contacts/rmsd/"+path+"/"+ "*rmsf.txt")
df1=pd.concat((pd.read_csv(infile1,sep='\t',header=None, names=['Pos',infile1.split("/")[-1].split("_")[0]+"_"+infile1.split("/")[-1].split("_")[1]]) for infile1 in all_files), axis=1)
df1 = df1.loc[:,~df1.columns.duplicated(keep='first')]
df1=df1.set_index('Pos')
df2=df1.T
path1='_'.join(path.split('_')[1:5])
df3=pd.read_csv("./plots/violin_Gs_Gi_allrmsd_"+path1+".tsv", comment="#", sep="\t")
df3=df3[['pdb1', 'pdb2','gfam1']].drop_duplicates()
df2=df2.reset_index()
df2[['pdb1', 'pdb2']] = df2['index'].str.split('_', 1, expand=True)
df4 = df3.merge(df2,left_on=['pdb1','pdb2'],right_on=['pdb1','pdb2'],how="left")
df4=df4.drop_duplicates()
df4=df4.drop(['index'], axis=1)
df4.to_csv(path1+'_rmsf.tsv',index=False,header=True,sep='\t')
    
    