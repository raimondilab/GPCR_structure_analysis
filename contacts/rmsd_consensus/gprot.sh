#!/bin/bash
awk -F '\t' '{printf("%s_%s\n",$1,$9)}' GPCR_structs_clean_consensus.tsv > Gprot_chains.txt
sed -i '1d' Gprot_chains.txt
gunzip ./cifs/*.gz
ls ./cifs/*.cif | xargs -n 1 ./test.py
./split_fasta.py < Gprot_chains.txt > ./seqs/Gprot.fa
clustalo  --hmm-in=G-alpha.hmm  -i ./seqs/Gprot.fa -o ./seqs/Gprot_ali.fa --outfmt=fasta
./msa_consensus_mappdb.py ./seqs/Gprot_ali.fa > Gprot_pdbs_consensus.txt 
