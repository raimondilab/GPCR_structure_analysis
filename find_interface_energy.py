#!/usr/bin/env python3

import csv,sys

# Parse Rosetta output files and appends the interface energy to the list of structures

filein=open(sys.argv[1])
read_tsv=csv.reader(filein, delimiter="\t")
fout=open(sys.argv[2],'wt')
write_tsv=csv.writer(fout, delimiter="\t")
f=0
for row in read_tsv:
	if f==0:
		f+=1
		write_tsv.writerow(row+['Interface_energy'])
		continue
	try:
		fileen=open("../GPCR_experimental_structures/relaxed_struct/"+row[0]+'_interface.tsv')
	except:
		print(row[0])
		continue
	read_en=csv.reader(fileen, delimiter=" ")
	flag=0
	for line in read_en:
		flag+=1
		if flag==3:
			k=0
			h=-1
			while k<6: #We take the 6th field: 
				h+=1
				if line[h]!='':
					k+=1
			write_tsv.writerow(row+[line[h]])
			break
	fileen.close()
filein.close()
fout.close()
