#!/usr/bin/env python3
import pandas as pd
import os,sys,operator, math
import gzip, sqlite3
import numpy as np
import glob, itertools
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import *
from difflib import SequenceMatcher

### extract fasta from pdbs
df=pd.read_csv("./use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
pdbs=[]
for index,row in df.iterrows():
	pdbs.append(row['PDB_ID']+"_"+row['Gprotein_chain'])


####I use here the same workaround employed for XL_interpreter when applied to user input structures (i.e. not taken from the PDB)
### I create a db for blast searches with the sequences from the input PDB and I blast on hit SPROT sequences 
#pdbdir="/net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz" % (pdb[1:3],pdb)
for i in pdbs:
    pdbdir='./cifs/'+ i.split('_')[0]+'.cif.gz'
    pdbpath=gzip.open(pdbdir, "rt")
    print(pdbpath)
    seq_list=[]
    for record in SeqIO.parse(pdbpath, "cif-atom"):
        pdb=i.split('_')[0].strip()
        print (pdb+"_"+record.annotations["chain"], record.seq)
        seq=SeqRecord(record.seq, id=pdb+"_"+record.annotations["chain"], description="Sequence from PDB")
        seq_list.append(seq)
        SeqIO.write(seq_list,"./seqs/"+"%s.fa" % (pdb), "fasta")

### combine gprotein fasta       

with open('./seqs/Gprot.fa','w') as f:
    for infile in glob.glob("./seqs/*.fa"):
    #print (infile)
        infasta=open(infile, "r")
        infasta2=infasta.read()
        saved=""
        for fa in infasta2.split(">")[1:]:
            ID = fa.split(" ")[0]
            if ID in pdbs:
                saved=">"+fa
                print (saved,file=f)
      
 ### do alignment
command='clustalo  --hmm-in=G-alpha.hmm  -i ./seqs/Gprot.fa -o ./seqs/Gprot_ali.fa --outfmt=fasta'
os.system(command)

fs="\t"


###Script to generate consensus columns from an input MSA generated from PDB sequences
### and map back MSA consensus positions to their counterparts in the PDB residue numbers (i.e. _atom_site.label_seq_id obtained through the biopython residue method .get_id())


track={}
aln = AlignIO.read('./seqs/Gprot_ali.fa', "fasta")
ref_flag=0
    
    
#print ("REFSEQ", ref_id)
#print ("Seq_id\tDeletions\tInsertions\tDivergent_positions")
gaps=[".", "-", "X"]
poscc={}
for a in aln:
  aid=a.id.split(",")[0]
  track[aid]={}
  ref_pos=0
  tar_pos=0 
  ii=0
  jj=0
  for aa in a.seq:
    if aa not in gaps:
      track[aid][ii]=jj
      jj+=1
      if ii not in poscc:
        poscc[ii]=1
      else:
        poscc[ii]+=1
    ii+=1


conpos=[]
for k in poscc.keys():
  if poscc[k] == len(aln):
    conpos.append(k)

with open('Gprot_pdb_consensus.txt', 'w') as outfile:
    for aid in track.keys():
        outstr = aid + fs
        pdb = aid.split("_")[0]
        chain = aid.split("_")[1]
        pdbdir = "./cifs/" + pdb + ".cif.gz"
        pdbpath = gzip.open(pdbdir, "rt")
        parser = PDB.MMCIFParser(QUIET=True)
        struct = parser.get_structure('XXX', pdbpath)
        model = struct[0]
        pdbchain = {}
        for ii, res in enumerate(model[chain]):
            pdbchain[ii] = res.get_id()[1]
        for msap in conpos:
            p = track[aid][msap]
            outstr += str(pdbchain[p]) + ","
        outfile.write(outstr + '\n')
