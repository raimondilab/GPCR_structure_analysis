#!/data/SW/anaconda3/envs/myenv/bin/python

import numpy as np
import sys
import os

input_file=sys.argv[1]
pdb=input_file.split("/")[2].split(".")[0]

ids={}
with open('./use_file/GPCR_structs_clean_consensus.tsv') as e:
    for line in e:
        pdb1=line.split('\t')[0]
        uniprot=line.split('\t')[1]
        gpcr = line.split('\t')[2]
        gpcr_chain=line.split('\t')[4]
        gprot=line.split('\t')[6]
        gprotid=line.split('\t')[7]
        gprot_chain=line.split('\t')[8].strip('\n')
        if pdb in pdb1:
            ids[pdb]=[gpcr,uniprot,gprot,gprotid,gpcr_chain,gprot_chain]
            output='./cont_file/'+pdb+"_"+gpcr+"_"+gprot+"_cons1.txt"

with open(output, 'w') as f1:
    print('GPCR'+'\t'+'Uniprot'+'\t'+'Pos1'+'\t'+'Gprotein'+'\t'+'Gprotein_id' + '\t'+'Pos2',file=f1)
    with open(input_file) as f:
        for line in f:
            chain1 = line.split(' ')[0].split('/')[0]
            chain2=line.split(' ')[1].split('/')[0]
            num1=line.split(' ')[0].split('/')[1]
            num2=line.split(' ')[1].split('/')[1]
            for k,v in ids.items():
                 if v[4] in chain1 and v[5] in chain2:
                    print(v[0]+'\t'+v[1]+'\t'+num1+'\t'+v[2] +'\t'+v[3]+'\t'+num2,file=f1)
                 elif  v[4] in chain2 and v[5] in chain1:
                     print(v[0]+'\t'+v[1]+'\t'+num2+'\t'+v[2] +'\t'+v[3]+'\t'+num1,file=f1)
