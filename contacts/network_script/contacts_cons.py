#!/data/SW/anaconda3/envs/myenv/bin/python

import sys
import os
import sqlite3

input_file=sys.argv[1]
pdb=input_file.split("/")[2].split(".")[0]
ids = []
with open('./use_file/GPCR_structs_clean_consensus.tsv') as e:
    for line in e:
        l = line.split('\t')
        pdb1=l[0]
        gpcr_uniprot=l[1]
        gpcr = l[2]
        gpcr_chain=l[4]
        gprot_uniprot = l[5]
        gprot=l[6]
        gprotid = l[7]
        gprot_chain=l[8].strip('\n')
        if pdb == pdb1:
            ids = [gpcr, gpcr_uniprot, gprot, gprotid, gpcr_chain, gprot_chain, gprot_uniprot]
            output="./cont_file/"+pdb+"_"+gpcr+"_"+gprot+"_cont.txt"
            break
            
# Create the pdb to uniprot mapping dictionary and check if pdb has been found in GPCR_structs_clean.tsv
try:
    pdb_uniprot = {ids[4]:{}, ids[5]:{}}
except IndexError:
    print(pdb+" is not listed in GPCR_structs_clean.tsv")
    sys.exit(0)

# Find the mapping using the database generated form SIFTS
conn = sqlite3.connect("/data/DB/SIFTS/pdb2uniprot_mappings.db")
c = conn.cursor()
for line in c.execute('select * from pdb2uniprot where pdbid = '+"'"+ pdb+"'"):
    if line[2] in [ids[1], ids[6]]:  # Check if the uniprot accession is the right one
        pdbpos = line[1].split('|')
        if pdbpos[0] in ids[4:6]:  # Check if the chain is the right one
            pdb_uniprot[pdbpos[0]][pdbpos[1][3:]] = line[3][1:]

unmapped = 0
with open(output, 'w') as f1:
    print('GPCR'+'\t'+'Uniprot'+'\t'+'Pos1'+'\t'+'Gprotein'+'\t'+'Gprotein_id' + '\t'+'Pos2',file=f1)
    with open(input_file) as f:
        for line in f:
            chain1 = line.split(' ')[0].split('/')[0]
            chain2=line.split(' ')[1].split('/')[0]
            num1=line.split(' ')[0].split('/')[1]
            num2=line.split(' ')[1].split('/')[1]
            if ids[4] == chain1 and ids[5] == chain2:
                try:  # Conversion, from pdb to uniprot number, it might fail if it's an unmapped position 
                    num1 = pdb_uniprot[chain1][num1]
                    num2 = pdb_uniprot[chain2][num2]
                except KeyError:  # Avoids printing unmapped positions, but keeps track of the problem
                    unmapped += 1
                    continue
                print(ids[0]+'\t'+ids[1]+'\t'+num1+'\t'+ids[2] +'\t'+ids[3]+'\t'+num2,file=f1)
            elif ids[4] == chain2 and ids[5] == chain1:
                try:  # Same as a few lines before
                    num1 = pdb_uniprot[chain1][num1]
                    num2 = pdb_uniprot[chain2][num2]
                except KeyError:
                    unmapped += 1
                    continue
                print(ids[0]+'\t'+ids[1]+'\t'+num2+'\t'+ids[2] +'\t'+ids[3]+'\t'+num1,file=f1)

if unmapped > 0:  # Reports the number of unmatched contacts if present
    print(pdb+ " has "+ str(unmapped)+" unmapped contacts.")