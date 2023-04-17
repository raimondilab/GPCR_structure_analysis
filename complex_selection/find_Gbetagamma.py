#!/usr/bin/python3

'''This module adds information about the presence of Gbeta-gamma in the structure'''

import gzip
import csv
import sqlite3
import sys
import pandas as pd

conn = sqlite3.connect("/projects/bioinformatics/DB/uniprot/Uniac2GN.db")
c = conn.cursor()

# Read the dataframe sys.argv[1]
GPCR_structs = pd.read_csv(sys.argv[1], sep='\t')
structs=set(GPCR_structs['PDB_ID'])

G_beta_gamma = {}
for pdb in structs:
    G_beta_gamma[pdb] = ["-"]*8

# We use the InterPro entries IPR016346 and IPR036284 respectively for Gbeta and Ggamma
G_beta_names = {}
G_gamma_names = {}
fin = gzip.open('/projects/bioinformatics/DB/SIFTS/pdb_interpro.tsv.gz', 'rt')
read_tsv = csv.reader(fin, delimiter="\t")
for row in read_tsv:
    if row[0] in structs:
        if row[3].split('_')[0] == 'IPR016346':
            if row[2] not in G_beta_names:
                for line in c.execute("SELECT * FROM uniac2gn WHERE uniac = '"+row[2]+"'"):
                    G_beta_names[row[2]] = [line[2].upper(), line[1]]
            G_beta_gamma[row[0]][:4] = [row[1]]+[row[2]]+G_beta_names[row[2]]
        if row[3].split('_')[0] == 'IPR036284':
            if row[2] not in G_gamma_names:
                for line in c.execute("SELECT * FROM uniac2gn WHERE uniac = '"+row[2]+"'"):
                    G_gamma_names[row[2]] = [line[2].upper(), line[1]]
            G_beta_gamma[row[0]][4:] = [row[1]]+[row[2]]+G_gamma_names[row[2]]
fin.close()
conn.close()

# Create new columns in GPCR_structs containing the information about Gbeta-gamma
GPCR_structs['Gbeta_chain'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][0])
GPCR_structs['Gbeta_Uniprot_AC'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][1])
GPCR_structs['Gbeta_Gene_name'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][2])
GPCR_structs['Gbeta_Uniprot_ID'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][3])
GPCR_structs['Ggamma_chain'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][4])
GPCR_structs['Ggamma_Uniprot_AC'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][5])
GPCR_structs['Ggamma_Gene_name'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][6])
GPCR_structs['Ggamma_Uniprot_ID'] = GPCR_structs['PDB_ID'].map(lambda x: G_beta_gamma[x][7])

# Write the new dataframe
GPCR_structs.to_csv(sys.argv[2], sep='\t', index=False)
