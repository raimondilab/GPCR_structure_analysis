#!/data/SW/anaconda3/envs/myenv/bin/python

import sys, os, operator,gzip,glob, math
fs='\t'
path=sys.argv[1]
for infile in glob.glob(path+"*fit.txt"):
    pdb1=infile.split("/")[-1].split("_")[0]
    pdb2=infile.split("/")[-1].split("_")[1]
    for l in open(infile, "r"):
        if  l.find("Selection RMSD") != -1:
            rms=float(l.split()[3].strip("\n"))
            print (pdb1,fs,pdb2,fs,rms)
