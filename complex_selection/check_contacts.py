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
                unmapped += 1
                continue
            write_cont.writerow([ids[2]]+[ids[1]]+[num1]+ids[6:8]+[num2])
        elif ids[4] == chain2 and ids[8] == chain1:
            try:  # Same as a few lines before
                num1 = pdb_uniprot[chain1][num1]
                num2 = pdb_uniprot[chain2][num2]
            except KeyError:
                unmapped += 1
                continue
            write_cont.writerow([ids[2]]+[ids[1]]+[num2]+ids[6:8]+[num1])
    fin.close()
    out.close()
    if unmapped > 0:  # Reports the number of unmatched contacts if present
        print(ids[0]+ " has "+ str(unmapped)+" unmapped contacts.")
        return False
    return True
    
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
promiscuous = {"RHO": "Yes",
                "CALCRL": "No",
                "GPR52": "No",
                "SMO": "Yes",
                "SCTR": "No",
                "GABBR2": "No",
                "GPR88": "No", 
                "FZD7": "Yes", 
                "PTH2R": "Yes", 
                "MC1R": "No", 
                "CCR1": "No", 
                "VIPR2": "No", 
                "GPR139": "Yes", 
                "ADGRG2": "Yes", 
                "NPY2R": "Yes", 
                "NPY4R": "Yes", 
                "CX3CR1": "No", 
                "MRGPRD": "Yes",
                "ADGRD1": "Yes",
                "ADGRG5": "No",
                "ADGRL3": "Yes",
                "ADGRG1": "Yes",
                "GHRHR": "No",
                "ADGRF1": "Yes",
                "ADGRG4": "No",
                "HTR5A": "Yes",
                "CCR3": "No",
                "CCR2": "No",
                "TSHR": "Yes",
                "MC2R": "No"}
filein = open("meta_encoded.txt")
filein.readline()
read_tsv = csv.reader(filein, delimiter="\t")
for row in read_tsv:
    if row[5].count('1') > 1:
        promiscuous[row[0]] = "Yes"
    else:
        promiscuous[row[0]] = "No"

# Keep only complexes for which all the contacts are mapped
filein = open("GPCR_structs.tsv")
read_tsv = csv.reader(filein, delimiter="\t")
fout = open("GPCR_structs_clean.tsv", "wt")
write_tsv = csv.writer(fout, delimiter='\t')
flag = 0
for row in read_tsv:
    if flag == 0:
        flag += 1
        write_tsv.writerow(row+["Antibody", "Promiscuous", "Class"])#, "Fingerprint"])
        continue
    if map_contact(row):
        write_tsv.writerow(row+[antibody(row[0]), promiscuous[row[2]], classification[row[2]]])#, promiscuity[row[2]]])
filein.close()
fout.close()
conn.close()
