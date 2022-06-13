#!/data/SW/anaconda3/envs/myenv/bin/python

import pandas as pd
import numpy as np
import sys
import os

input_file=sys.argv[1]
pdb=input_file.split("/")[2].split("_")[0]
gpcr1=input_file.split("/")[2].split("_")[1]
gprot2=input_file.split("/")[2].split("_")[2]
output="./ann_file/"+pdb+"_"+gpcr1+"_"+gprot2+"_2.txt"
gprot=pd.read_csv("./use_file/human_sequence_table_FR_new.txt", comment="#", sep="\t",index_col=0)
gprot1=gprot.unstack().reset_index()
gprot1.columns=['Gprotein_id','CGN','Pos2']
gprot1['Pos2'] = gprot1['Pos2'].str.replace(r'\D', '')
gprot1['Pos2'] .replace('', np.nan, inplace=True)
gprot1=gprot1.dropna()
gpcr=pd.read_csv("./use_file/output_gpcrdb.tsv", dtype=object, sep="\t")
df=pd.read_csv(input_file, sep="\t", dtype=object)
d2=df.merge(gprot1,left_on=['Gprotein_id','Pos2'],right_on=['Gprotein_id','Pos2'],how="left")
d3=d2.merge(gpcr,left_on=['Uniprot','Pos1'],right_on=['Uniprot','Pos1'],how="left")
d3.to_csv(output,sep='\t',index=False)
