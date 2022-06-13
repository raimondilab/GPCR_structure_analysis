#!/data/SW/anaconda3/envs/myenv/bin/python

from collections import defaultdict
import pandas as pd
import glob

df=pd.read_csv("./use_file/GPCR_structs_clean_consensus.tsv", comment="#", sep="\t")
s=df.groupby(['Receptor_Gene_name','Gprotein_Gene_name'])['PDB_ID'].apply(list).to_dict()
all_files = glob.glob( "./ann_file/*2.txt")    
fil={}
for i in all_files:
    pdb=i.split("/")[2].split("_")[0]
    gpcr=i.split("/")[2].split("_")[1]
    gprot=i.split("/")[2].split("_")[2]
    fil[i]=(gpcr,gprot)

swapped = defaultdict(set)
for k, v in fil.items():
    swapped[v].add(k)           

alls={}           
for k,v in swapped.items():
    dfs=[]
    for i in v:
      dfs.append(pd.read_csv(i,sep='\t'))
      alls[k]=pd.concat(dfs)
      alls[k]=alls[k].drop(alls[k].columns[[1,4,7]], axis=1).drop_duplicates().dropna()
     
        
for k,v in alls.items():
    v.to_csv('./gp_pairs/'+k[0]+"_"+k[1]+".tsv",index=None,header=True,sep='\t')

alld=pd.concat(alls.values(), ignore_index=True)
alld.to_csv("./use_file/all_gpcrs_cons.tsv",index=None,header=True,sep='\t')


