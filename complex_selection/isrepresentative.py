#!/usr/bin/python3

''' This module is used to add a column to the tsv file to indicate whether the
structure is among the consensus structures.'''

import sys
import csv

# Read the consensus structures
filein = open(sys.argv[1])
filein.readline()
read_tsv = csv.reader(filein, delimiter="\t")
consensus = set()
for row in read_tsv:
    consensus.add(row[0])
filein.close()

# Read the list of structures and add a column to indicate whether the structure is in the consensus set
filein = open(sys.argv[2])
flag = 0
read_tsv = csv.reader(filein, delimiter="\t")
fout = open(sys.argv[3], 'wt')
write_tsv = csv.writer(fout, delimiter='\t')
for row in read_tsv:
    if flag == 0:
        flag += 1
        write_tsv.writerow(row + ["Representative"])
        continue
    if row[0] in consensus:
        write_tsv.writerow(row + ["Yes"])
    else:
        write_tsv.writerow(row + ["No"])
filein.close()
fout.close()
