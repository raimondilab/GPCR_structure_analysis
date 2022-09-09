#!/data/SW/anaconda3/envs/myenv/bin/python

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
sys.path.insert(0, "/data/Users/francesco/code/")
from sequence_tools import ExtractFasta, GetBlastOut, GetBlastOutX
f="\t"


def GetFasta(infile, outseqs):
  seq=""
  flag=0
  tracker=[]
  if infile[:-3] == ".gz":
    infile_stream=gzip.open(infile,"r")
  else:
    infile_stream=open(infile, "r")
  for l in infile_stream:
    if l[0] == ">":
      flag=0
      seq=""
      #gn=l.split()[0]
      seqid=l.split()[0].lstrip(">")
      outseqs[seqid]=l
    if l[0] != ">":
      outseqs[seqid]+=l.strip("\n")
	  
	  
pdb=sys.argv[1]
pdbid=pdb.split("/")[2].split('.')[0]
pdbid=pdbid.lower()


###Getting uniprot accession corresponding to each PDB chain
unip_list=[]
pdb2uac={}
for l in gzip.open("/data/DB/SIFTS/pdb_chain_uniprot.tsv.gz", "rt"):
  if l[0] != "#":
    pid=l.split(f)[0]
    chain=l.split(f)[1]
    uac=l.split(f)[2]
    if pid == pdbid:
      if pid not in pdb2uac:
        pdb2uac[pid]=[[chain,uac]]
      elif [chain, uac] not in pdb2uac[pid]:
        pdb2uac[pid].append([chain, uac])
      unip_list.append(uac)


####Retrieving corresponding fasta sequences
fastaname=("uniprot_query.fasta")
pac2pid, pac2gn,origid2pac = ExtractFasta("/data/DB/uniprot/uniprot_seq_ids_new.db", unip_list, fastaname)


####I use here the same workaround employed for XL_interpreter when applied to user input structures (i.e. not taken from the PDB)
### I create a db for blast searches with the sequences from the input PDB and I blast on hit SPROT sequences 
#pdbdir="/net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz" % (pdb[1:3],pdb)
pdbdir="%s" % (pdb)
pdbpath=open(pdbdir, "rt")
print(pdbpath)
seq_list=[]
for record in SeqIO.parse(pdbpath, "cif-atom"):
  pdb=pdb.strip()
  print (pdbid+"_"+record.annotations["chain"], record.seq)
  seq=SeqRecord(record.seq, id=pdbid+"_"+record.annotations["chain"], description="Sequence from PDB")
  seq_list.append(seq)

SeqIO.write(seq_list,"./seqs/"+"%s.fa" % (pdbid), "fasta")
