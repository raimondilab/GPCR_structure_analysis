#!/usr/bin/python3

'''This module outputs a file with a list of GCPR-Galpha pairs appearing in the same
structure on the pdb, as a result some structures may appear multiple times'''

import gzip
import csv
import sqlite3

# Read the list of Galpha containing structures

conn = sqlite3.connect("/home/pmiglionico/Uniac2GN.db")
c = conn.cursor()

filein = open("Galpha_list.txt")
read_tsv = csv.reader(filein, delimiter="\t")
Galpha = {}
Galpha_names = {}
for row in read_tsv:
    if row[2] not in Galpha_names:
        for line in c.execute("SELECT * FROM uniac2gn WHERE uniac = '"+row[2]+"'"):
            Galpha_names[row[2]] = [line[2].upper(), line[1]]
    if row[0] not in Galpha:
        Galpha[row[0]] = []
    Galpha[row[0]].append([row[0]]+[row[2]]+Galpha_names[row[2]]+[row[1]])
filein.close()

# Find GPCR containing chains

filein = gzip.open("/home/pmiglionico/pdb_chain_pfam.tsv.gz", 'rt')
read_tsv = csv.reader(filein, delimiter="\t")
GPCR = {"PF00001", "PF00002", "PF00003", "PF01534"}  # all Pfam entries corresponding to GPCRs
GPCR_names = {}
flag = 0
pairs = []
virus = {"HHV8P", "EBVB9", "HCMVA", "HCMVT", "RCMVM", "HHV6U", "MUHVS", "HHV7J", "HCMVM", "HHV6Z"}
for row in read_tsv:
    if flag == 0: # skip first row
        flag = 1
        continue
    if row[3] in GPCR and row[0] in Galpha:
        if row[2] not in GPCR_names:
            for line in c.execute("SELECT * FROM uniac2gn WHERE uniac = '"+row[2]+"'"):
                GPCR_names[row[2]] = [line[2].upper(), line[1]]
        if GPCR_names[row[2]][1].split('_')[1] in virus:
            continue
        for el in Galpha[row[0]]:
            # Reorder the output entries to make them coherent with older files made by Marin
            pairs.append([el[0]]+[row[2]]+GPCR_names[row[2]]+[row[1]]+el[1:])
filein.close()
conn.close()

#Write the output
fout = open("GPCR_structs.tsv", "wt")
write_tsv = csv.writer(fout, delimiter='\t')
write_tsv.writerow(["PDB_ID",
                    "Receptor_Uniprot_AC",
                    "Receptor_Gene_name",
                    "Receptor_Uniprot_ID",
                    "Receptor_Chain",
                    "Gprotein_Uniprot_AC",
                    "Gprotein_Gene_name",
                    "Gprotein_Uniprot_ID",
                    "Gprotein_chain",
                    "Receptor_organism",
                    "Gprotein_organism"])
for el in pairs:
    write_tsv.writerow(el+[el[3].split('_')[1], el[7].split('_')[1]])
fout.close()
