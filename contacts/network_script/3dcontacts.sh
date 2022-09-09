#!/bin/bash


###The 1st taks is computing the residue contacts 
for pdb in `ls -1 ./cifs/*cif.gz`
do
  gunzip $pdb
  echo "${pdb%%.cif.gz}.cif"
  /data/Users/francesco/code/cpp/cifparse-obj-v7.105-prod-src/parser-test-app/bin/3DContact ${pdb%%.cif.gz}.cif 8 > $pdb\_contacts.txt
done

/data/Users/francesco/code/3Dcontact2Uniprot_via_SIFTsxml.py

gzip /data/Users/marin/contact_analysis/consensus/cifs/*cif
