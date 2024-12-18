#!/usr/bin/env python3

'''This module takes the computed iutput of the Interface energy
analysis and adds it to the list of GPCR-Galpha complexes'''

import csv
import sys

filein = open(sys.argv[1])
read_tsv = csv.reader(filein, delimiter="\t")
fout = open(sys.argv[2], 'wt')
write_tsv = csv.writer(fout, delimiter="\t")
f = 0
for row in read_tsv:
    if f == 0:
        f += 1
        write_tsv.writerow(row+['Interface_energy', 'dG/dSASA', 'dSASA'])
        continue
    try:
        fileen = open("../GPCR_experimental_structures/relaxed_struct_no_Ab/"+row[0]+'_interface.tsv')
    except FileNotFoundError:
        print(row[0])
        continue
    read_en = csv.reader(fileen, delimiter=" ")
    flag = 0
    for line in read_en:
        flag += 1
        if flag == 3:
            k = 0
            h = -1
            data = [0, 0]
            while k < 9:  #We take the 6th, the 7th and the 9th field
                h += 1
                if line[h] != '':
                    k += 1
                    if k == 6:
                        data[0] = line[h]
                    elif k == 7:
                        data[1] = line[h]
            write_tsv.writerow(row+data+[line[h]])
            break
    fileen.close()
filein.close()
fout.close()
