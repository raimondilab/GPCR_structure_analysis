#!/data/SW/anaconda3/envs/myenv/bin/python

'''This module further cleans the list generated bu clean_selection.py removing the complexes
in which not all contacts are mapped to the same Uniprot entry, this avoids chimeric complexes.
Moreover it generates a file for each PDB containing the mapping of the contacts to uniprot sequence'''

import csv
import sqlite3

# Connecting to the PDB-Uniprot residue-level mapping database
conn = sqlite3.connect("/home/pmiglionico/pdb2uniprot_mappings.db")
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
    out = open("../GPCR_experimental_structures/cont_file/"+ids[0]+"_"+ids[2]+"_"+ids[6]+"_cont.txt", "w")
    write_cont = csv.writer(out, delimiter='\t')
    write_cont.writerow(['GPCR', 'Uniprot', 'Pos1', 'Gprotein', 'Gprotein_id', 'Pos2'])
    fin = open("/home/pmiglionico/pdb-mmCIF_CBcontacts_blast/"+ids[0][1:3]+"/"+ids[0]+"_3dc.txt")
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
            write_cont.writerow(ids[0:2]+[num1]+ids[6:7]+[num2])
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

filein = open("GPCR_structs.tsv")
read_tsv = csv.reader(filein, delimiter="\t")
fout = open("GPCR_structs_clean.tsv", "wt")
write_tsv = csv.writer(fout, delimiter='\t')
flag = 0
# Keep only compleces for which all the contacts are mapped
for row in read_tsv:
    if flag == 0:
        flag += 1
        write_tsv.writerow(row)
        continue
    if map_contact(row):
        write_tsv.writerow(row)
filein.close()
fout.close()
conn.close()
