#!/bin/bash

zgrep PF00503 /home/fraimondi/BIOINFO1_DB/SIFTS/pdb_chain_pfam.tsv.gz > Galpha_list.txt
python3 find_GPCR.py
rm Galpha_list.txt
