#!/usr/bin/env python3


import os,sys,operator
import gzip, re, glob
import subprocess

def GetFasta(infile, idlist, outseqs):
  seq=""
  flag=0
  for l in gzip.open(infile, "rt"):
    if l[0] == ">":
      flag=0
      if seq != "":
        outseqs.append(seq)
      seq=""
      uniac=l.split("|")[1]
      for chunk in l.split():
        if chunk.find("GN=") != -1:
          gn=chunk.split("=")[1]
      if l.find("HUMAN") != -1:
        if uniac in idlist:
          flag=1
        elif gn in idlist:
          flag=1
        else:
          continue
    if flag == 1:
      seq+=l
    else:
      continue
        
        
    

id_list=[]
for l in sys.stdin:
  id_list.append(l.strip("\n"))
  
for infile in glob.glob("./seqs/*.fa"):
  #print (infile)
  infasta=open(infile, "r")
  infasta2=infasta.read()
  saved=""
  for fa in infasta2.split(">")[1:]:
    ID = fa.split("\n")[0].split()[0]
    if ID in id_list:
      saved=">"+fa
  print (saved)
      
      
