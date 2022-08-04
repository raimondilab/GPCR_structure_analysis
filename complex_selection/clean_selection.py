#!/usr/bin/python3

'''This module takes a list of GPCR-Galpha pairs produced by find_GPCR.py and elimminates
redundancies, considering only the most representative GPCR-Galpha pair for each structure'''

import csv
import sqlite3
import sys

# Read the list of structures
filein = open("GPCR_structs.tsv")
read_tsv = csv.reader(filein, delimiter="\t")
pair = {}
flag = 0
for row in read_tsv:
    if flag == 0:
        flag += 1
        header = row
        continue
    if row[0] not in pair:
        pair[row[0]] = []
    pair[row[0]].append(row)
filein.close()

# Find the resolution of each structure
resolution = {}
for pdb in pair:
    filein = open("../GPCR_experimental_structures/structures/"+pdb+".cif")
    for row in filein:
        if row.find("resolution") != -1:
            h = row.split(' ')
            if (h[0][-len("resolution"):] == "resolution" or
                    h[0][-len("resolution_high"):] == "resolution_high"):
                k = 1
                while h[k] == '':
                    k += 1
                resolution[pdb] = float(h[k])
                break
    filein.close()

# Retrieve the number of contacts between each chain pair
for pdb in pair:
    try:
        filein = open("/home/pmiglionico/pdb-mmCIF_CBcontacts_blast/"+pdb[1:3]+"/"+pdb+"_3dc.txt")
    except FileNotFoundError:
        sys.exit("Error: please update contacts,"+pdb+" not found")
    read_cont = csv.reader(filein, delimiter=" ")
    contact = {}
    for row in read_cont:
        if row == []:
            break
        ch1 = row[0].split('/')[0]
        ch2 = row[1].split('/')[0]
        if ch1 == ch2:
            continue
        if ch1 not in contact:
            contact[ch1] = {}
        if ch2 not in contact[ch1]:
            contact[ch1][ch2] = 0
        contact[ch1][ch2] += 1
    filein.close()
    for el in pair[pdb]:
        try:
            el.append(contact[el[4]][el[8]])
        except KeyError:
            try:
                el.append(contact[el[8]][el[4]])
            except KeyError:
                el.append(0)

# Find the coverage of each chain
conn = sqlite3.connect("/home/pmiglionico/pdb2uniprot_mappings.db")
c = conn.cursor()
for pdb in pair:
    mapping = []
    for row in c.execute('select * from pdb2uniprot where pdbid = '+"'"+pdb+"'"):
        mapping.append(row[1:3])
    for el in pair[pdb]:
        el += [0, 0]
        for res in mapping:
            ch = res[0].split('|')[0]
            if ch == el[4] and res[1] == el[1]:
                el[-2] += 1
            if ch == el[8] and res[1] == el[5]:
                el[-1] += 1
conn.close()

# Find the length of each chain fo find the coverage ratio
conn = sqlite3.connect("/home/fraimondi/BIOINFO1_DB/uniprot/uniprot_seq_ids_new.db")
c = conn.cursor()
for pdb in pair:
    for el in pair[pdb]:
        el += [0, 0, resolution[pdb]]
        for line in c.execute('select * from uniprot where id_ac = '+"'"+ el[1]+"'"):
            el[-3] = el[-5]/len(''.join(line[3].split('\n')[1:]))
        for line in c.execute('select * from uniprot where id_ac = '+"'"+ el[5]+"'"):
            el[-2] = el[-4]/len(''.join(line[3].split('\n')[1:]))
conn.close()

# Write the output
header += ["Contacts",
           "Receptor_residues",
           "Gprotein_residues",
           "Receptor_coverage",
           "Gprotein_coverage",
           "Resolution"]
fout = open("GPCR_structs_clean.tsv", "wt")
write_tsv = csv.writer(fout, delimiter='\t')
write_tsv.writerow(header)
for pdb in pair:
    # Select the pair of chains with the highest number of contacts (to solve the problem of dimers
    # The pair of Uniprot accessions with the highest coverage is used to name the two chains
    pair[pdb].sort(key=lambda x: (x[-6], x[-3]*x[-2]))
    write_tsv.writerow(pair[pdb][-1])
fout.close()
