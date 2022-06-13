#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys,operator, math, glob, pickle
import gzip
import numpy as np
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from scipy.stats import hypergeom
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import *



fs="\t"


###Script to generate consensus columns from an input MSA generated from PDB sequences
### and map back MSA consensus positions to their counterparts in the PDB residue numbers (i.e. _atom_site.label_seq_id obtained through the biopython residue method .get_id())

#print (ac2gn.keys())
###Loading the MSA through Biopython AlignIO function
track={}
aln = AlignIO.read(sys.argv[1], "fasta")
ref_flag=0
    
    
#print ("REFSEQ", ref_id)
#print ("Seq_id\tDeletions\tInsertions\tDivergent_positions")
gaps=[".", "-", "X"]
poscc={}
for a in aln:
  aid=a.id.split(",")[0]
  track[aid]={}
  ###Looping now through the sequences to spot differences
  ###Using three counters, one for the canonical (ref)sequence and one for the target (tar) sequence and one for the MSA
  ref_pos=0
  tar_pos=0 
  ii=0
  jj=0
  #print (aid, a.seq)
  for aa in a.seq:
    if aa not in gaps:
      track[aid][ii]=jj
      jj+=1
      #print (aa,jj)
      if ii not in poscc:
        poscc[ii]=1
      else:
        poscc[ii]+=1
    ii+=1


conpos=[]
for k in poscc.keys():
  if poscc[k] == len(aln):
    conpos.append(k)
#print (conpos)

  
for aid in track.keys():
  outstr=aid+fs
  pdb=aid.split("_")[0]
  chain=aid.split("_")[1]
  pdbdir="./cifs/"+pdb+".cif"
  pdbpath=open(pdbdir, "rt")
  parser = PDB.MMCIFParser(QUIET=True)
  struct = parser.get_structure('XXX',pdbpath)
  model = struct[0]
  pdbchain={}
  for ii, res in enumerate(model[chain]):
    #print (ii, res.get_id())
    pdbchain[ii]=res.get_id()[1]
  #print (aid, pdbchain)
  #print (track[aid])
  for msap in conpos:
    p=track[aid][msap]
    #if p in pdbchain:
    #print (aid, chain, msap, p)
    outstr+=str(pdbchain[p])+","
    #outstr+=str(p)+","
  print (outstr)
    
    


