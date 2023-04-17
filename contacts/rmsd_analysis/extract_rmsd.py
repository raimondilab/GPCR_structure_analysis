import sys
import os
import operator
import gzip
import glob
import math
#####Extract rmsd from the calculated files
fs = '\t'
path = sys.argv[1]
all_files = glob.glob(path + "/**/*_fit.txt", recursive=True)
for infile in all_files:
    pdb1 = infile.split("/")[-1].split("_")[0]
    pdb2 = infile.split("/")[-1].split("_")[1] 
    for l in open(infile, "r"):
        if l.find("Selection RMSD") != -1:
            rms = float(l.split()[3].strip("\n"))
            print(pdb1, fs, pdb2, fs, rms)

