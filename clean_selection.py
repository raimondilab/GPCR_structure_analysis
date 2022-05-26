#!/usr/bin/python3

import gzip,csv,sqlite3
from Bio.PDB import MMCIFParser

# Read the list of structures
filein=open("GPCR_structs.tsv")
read_tsv=csv.reader(filein, delimiter="\t")
pair={}
header=[]
flag=0
for row in read_tsv:
	if flag==0:
		flag+=1
		header=row
		continue
	if row[0] not in pair:
		pair[row[0]]=[]
	pair[row[0]].append(row)
filein.close()

# Find the resolution of each structure
resolution={}
mmCIF_parser = MMCIFParser(QUIET=True)
for pdb in pair:
	structure_cif = mmCIF_parser.get_structure(pdb,"../GPCR_experimental_structures/structures/"+pdb+".cif")
	resolution[pdb]=structure_cif.header["resolution"]

# Retrieve the number of contacts between each chain pair
for pdb in pair:
	filein=open("/home/pmiglionico/pdb-mmCIF_CBcontacts_blast/"+pdb[1:3]+"/"+pdb+"_3dc.txt")
	read_cont=csv.reader(filein, delimiter=" ")
	contact={}
	for row in read_cont:
		if row==[]:
			break
		ch1=row[0].split('/')[0]
		ch2=row[1].split('/')[0]
		if ch1==ch2:
			continue
		if ch1 not in contact:
			contact[ch1]={}
		if ch2 not in contact[ch1]:
			contact[ch1][ch2]=0
		contact[ch1][ch2]+=1
	for el in pair[pdb]:
		try:
			el.append(contact[el[4]][el[8]])
		except KeyError:
			el.append(contact[el[8]][el[4]])
		
# Find the coverage of each chain

conn = sqlite3.connect("/home/fraimondi/BIOINFO1_DB/SIFTS/pdb2uniprot_mappings.db")
c = conn.cursor()
for pdb in pair:
	mapping=[]
	for row in c.execute('select * from pdb2uniprot where pdbid = '+"'"+pdb+"'"):
		mapping.append(row[1:3])
	for el in pair[pdb]:
		el+=[0,0,resolution[pdb]]
		for res in mapping:
			ch=res[0].split('|')[0]
			if ch==el[4] and res[1]==el[1]:
				el[-3]+=1
			if ch==el[8] and res[1]==el[5]:
				el[-2]+=1
	
# Write the output
header+=["Contacts","Receptor_residues","Gprotein_residues","Resolution"]
fout=open("GPCR_structs_clean.tsv","wt")
write_tsv=csv.writer(fout, delimiter='\t')
write_tsv.writerow(header)
for pdb in pair:
	# Select the pair of chains with the highest number of contacts (to solve the problem of dimers
	# The pair of Uniprot accessions with the highest coverage is used to name the two chains
	pair[pdb].sort(key=lambda x:(x[-4],x[-3]*x[-2]))
	write_tsv.writerow(pair[pdb][-1])
fout.close()
