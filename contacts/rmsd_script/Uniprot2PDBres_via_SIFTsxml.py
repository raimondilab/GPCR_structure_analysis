#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys,operator, math,time
import gzip, sqlite3
import numpy as np
import glob, itertools
#from Bio import PDB
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.SeqRecord import SeqRecord
#from Bio import AlignIO
#from Bio.Align import AlignInfo
#from Bio.PDB import PDBIO
#from Bio.PDB.Polypeptide import *
#from difflib import SequenceMatcher
from Bio.SeqUtils import seq1
sys.path.insert(0, "/data/Users/francesco/code/")
#from sequence_tools import ExtractFasta
import xml.etree.ElementTree as ET
from os.path import exists

AA=['K','R','H','E','D','N','Q','Y','S','T','C','A','G','P','F','I','L','V','M','W']

fs="\t"

GAPS=[".", "-", "X"]


def UnfoldSiftXml(inpxml):
  mappings={}
  xml=gzip.open(inpxml, "rt")
  tree = ET.parse(xml)
  root = tree.getroot()
  for entity in root:
    #print(entity.tag, entity.attrib)
    for chain in entity:
      #print (chain.tag, chain.attrib)
      for listres in chain:
        #print (listres.tag, listres.attrib)
        for res in listres:
          #print (res.tag, res.attrib)
          flag1=0
          flag2=0
          for r in res:
            if r.attrib['dbSource'] == 'UniProt' and 'dbAccessionId' in r.attrib and 'dbResNum' in r.attrib and 'dbResName' in r.attrib:
              uniacc=r.attrib['dbAccessionId']
              uniresnum=r.attrib['dbResNum']
              uniresname=r.attrib['dbResName']
              flag1=1
            if r.attrib['dbSource'] == 'PDB' and 'dbAccessionId' in r.attrib and 'dbResNum' in r.attrib and 'dbResName' in r.attrib and 'dbChainId' in r.attrib:
              pdbacc=r.attrib['dbAccessionId']
              pdbresnum=r.attrib['dbResNum']
              pdbresname=r.attrib['dbResName']
              pdbchain=r.attrib['dbChainId']
              flag2=1
            if flag1 == 1 and flag2 == 1 and pdbresnum.find('null') == -1:
              if uniacc not in mappings:
                #mappings[pdbchain]={}
                #mappings[pdbchain][pdbresnum]=[pdbresname,uniacc,uniresnum,uniresname]
                mappings[uniacc]={}
                mappings[uniacc][uniresnum]=[pdbresname,pdbchain,pdbresnum,uniresname]
              elif uniacc in mappings:
                #mappings[pdbchain][pdbresnum]=[pdbresname,uniacc,uniresnum,uniresname]
                mappings[uniacc][uniresnum]=[pdbresname,pdbchain,pdbresnum,uniresname]
  return mappings




  

###Getting residue level mappings to Uniprot from the SIFT XML for a given PDB 
#for contact_file in glob.glob("*.gz_contacts.txt"):
pdbid=sys.argv[1]
xmlpath="/data/DB/SIFTS/ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/%s/%s.xml.gz" % (pdbid[1:3], pdbid)
if os.path.exists(xmlpath):
  uni2pdb=UnfoldSiftXml(xmlpath)
  for ac in uni2pdb.keys():
    for res in uni2pdb[ac].keys():
      chain=uni2pdb[ac][res][1]
      pdbresnum=uni2pdb[ac][res][2]
      print (pdbid+fs+chain+fs+pdbresnum+fs+ac+fs+res)
