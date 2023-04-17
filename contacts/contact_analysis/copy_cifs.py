import shutil
import pandas as pd
import os
##### Copy cifs files ,extracted contact interfaces for each cif file
df=pd.read_csv("./use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
pdb=df.PDB_ID.tolist()
names={}
for index,row in df.iterrows():
	names[row['PDB_ID']]=row['Receptor_Gene_name']+"_"+row['Gprotein_Gene_name']


pdb_pos=[]
for i in pdb:
    pdb_pos.append("/projects/bioinformatics/DB/pdb-mmCIF/"+i[1:3]+"/"+i+".cif.gz")


for j in pdb_pos:   
    dst_path = './cifs/'
    shutil.copy(j, dst_path)

for l in pdb:
    dst_path1= './cont_file/'
    if l =='5g53':
        print('No need')
    else:        
        shutil.copy("/home/pmiglionico/GPCR_experimental_structures/cont_file/"+l+"_cont.tsv", dst_path1)
     
 

input_files = [f for f in os.listdir('./cont_file/') if f.endswith('.tsv')]

for file in input_files:
    for k,v in names.items():
        if file.split('_')[0] == k:
            newName = "./cont_file/"+k+"_"+ v +"_cont.tsv"
            os.rename('./cont_file/'+file, newName)
           
 