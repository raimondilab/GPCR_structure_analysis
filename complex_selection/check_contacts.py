#!/usr/bin/python3

'''This module further cleans the list generated bu clean_selection.py removing the complexes
in which not all contacts are mapped to the same Uniprot entry, this avoids chimeric complexes.
Moreover it generates a file for each PDB containing the mapping of the contacts to uniprot sequence.
It also checks if an antibody is present in the structure and reports it in a new column in the table'''

import csv
import sqlite3
import Bio.PDB

# Connecting to the PDB-Uniprot residue-level mapping database
conn = sqlite3.connect("/projects/bioinformatics/DB/SIFTS/pdb2uniprot_mappings.db")
c = conn.cursor()

def map_contact(ids):

    '''This function checks the contacts for a single structure and generates the contact mapping file'''

    # Find the mapping using the database generated form SIFTS
    pdb_uniprot = {ids[4]:{}, ids[8]:{}}
    for line in c.execute('select * from pdb2uniprot where pdbid = '+"'"+ids[0]+"'"):
        if line[2] in [ids[1], ids[5]]:  # Check if the uniprot accession is the right one
            pdbpos = line[1].split('|')
            if pdbpos[0] in [ids[4], ids[8]]:  # Check if the chain is the right one
                pdb_uniprot[pdbpos[0]][pdbpos[1][3:]] = line[3][1:]

    flag = 0
    if ids[1] == "A0A0A0MTJ0": # Remap TSHR to the isoform present in SwissProt (there is only 1 mismatch)
        ids[1] = "P16473"
        ids[3] = "TSHR_HUMAN"
    if ids[1] == "Q9HB45": # Remap GHRHR to the isoform present in SwissProt (there is a 64 positions shift)
        ids[1] = "Q02643"
        ids[3] = "GHRHR_HUMAN"
        flag = 1

    # Check if ../GPCR_experimental_structures/cont_file/"+ids[0]+"_cont.tsv already exists
    # Since we corrected some of them manually we want to avoid overwriting them if they exist
    try:
        fin = open("../GPCR_experimental_structures/cont_file/"+ids[0]+"_cont.tsv")
        fin.close()
        return ids
    except FileNotFoundError:
        print("Mapping contacts for "+ids[0])
        pass

    # Map the contacts to Uniprot positions
    unmapped = 0
    out = open("../GPCR_experimental_structures/cont_file/"+ids[0]+"_cont.tsv", "w")
    write_cont = csv.writer(out, delimiter='\t')
    write_cont.writerow(['GPCR', 'Uniprot', 'Pos1', 'Gprotein', 'Gprotein_id', 'Pos2'])
    fin = open("/projects/bioinformatics/DB/pdb-mmCIF_CBcontacts/"+ids[0][1:3]+"/"+ids[0]+"_3dc.txt")
    for line in fin:
        if line == '\n':
            break
        chain1 = line.split(' ')[0].split('/')[0]
        chain2 = line.split(' ')[1].split('/')[0]
        num1 = line.split(' ')[0].split('/')[1]
        num2 = line.split(' ')[1].split('/')[1]
        if ids[4] == chain1 and ids[8] == chain2:
            try:  # Conversion from pdb to uniprot position, it fails if it's an unmapped position
                num1 = pdb_uniprot[chain1][num1]
                num2 = pdb_uniprot[chain2][num2]
            except KeyError:  # Avoids printing unmapped positions, but keeps track of the problem
                print(line.strip('\n'))
                unmapped += 1
                continue
            if flag == 1:
                num1 = int(num1) + 64
            write_cont.writerow([ids[2]]+[ids[1]]+[num1]+ids[6:8]+[num2])
        elif ids[4] == chain2 and ids[8] == chain1:
            try:  # Same as a few lines before
                num1 = pdb_uniprot[chain1][num1]
                num2 = pdb_uniprot[chain2][num2]
            except KeyError:
                print(line.strip('\n'))
                unmapped += 1
                continue
            if flag == 1:
                num2 = int(num2) + 64
            write_cont.writerow([ids[2]]+[ids[1]]+[num2]+ids[6:8]+[num1])
    fin.close()
    out.close()
    if unmapped > 0:  # Reports the number of unmatched contacts if present
        print(ids[0]+ " has "+ str(unmapped)+" unmapped contacts.")
    return ids
    
def antibody(pdb):
    
    '''This function finds nanobodies in the cif files'''
    
    structdict = Bio.PDB.MMCIF2Dict.MMCIF2Dict("../GPCR_experimental_structures/structures/"+pdb+".cif")
    for i in range(len(structdict['_entity.pdbx_description'])):
        descr = str.lower(structdict['_entity.pdbx_description'][i])
        # This list of if statements checks if in the entity desctiption there is the word "Antibody",
        # "Nanobody", "Fab" or "scF", it also checks for some typos that can be found here and there
        if descr.find('body') != -1:
            return descr
        elif descr.find('scf') != -1:
            return descr
        elif descr.find('fab') != -1:
            return descr
        elif descr.find('svf') != -1:
            return descr
        elif descr.find('nb') != -1:
            return descr
        elif descr.find('boy') != -1:
            return descr
    return "No"

classification = {}
filein = open("GPCRclassification.tsv")
filein.readline()
read_tsv = csv.reader(filein, delimiter="\t")
for row in read_tsv:
    classification[row[0]] = row[1]
filein.close()

# Read the list of PDBs to exclude, we compiled it manually
# The list contains PDBs that have problems in the mapping to uniprot sequences at the interface
filein = open("excluded_pdb.txt")
excluded = set()
for line in filein:
    excluded.add(line.strip('\n'))
filein.close()

# Keep only complexes for which all the contacts are mapped
filein = open("GPCR_structs.tsv")
read_tsv = csv.reader(filein, delimiter="\t")
fout = open("GPCR_structs_clean.tsv", "wt")
write_tsv = csv.writer(fout, delimiter='\t')
flag = 0
for row in read_tsv:
    if flag == 0:
        flag += 1
        write_tsv.writerow(row+["Antibody", "Class"])
        continue
    if row[0] not in excluded:
        row = map_contact(row)
        write_tsv.writerow(row+[antibody(row[0]), classification[row[2]]])
filein.close()
fout.close()
conn.close()
