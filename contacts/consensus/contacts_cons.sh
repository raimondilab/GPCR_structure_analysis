#!/bin/bash
./consensus_struct.py
awk -F '\t' '{print $1}' ./use_file/GPCR_structs_clean_consensus.tsv> pdb_list1.txt
for pdb in `cat pdb_list1.txt`; do echo "/data/DB/pdb-mmCIF/"${pdb:1:2}"/"$pdb".cif.gz" >> pdb_list_new1.txt ; done 
sed -i '1d' pdb_list_new1.txt
xargs -a pdb_list_new1.txt cp -t ./cifs/
./3dcontacts.sh 
ls ./cifs/*contacts.txt | xargs -n 1  ./contacts_cons.py
ls ./cont_file/*cons1.txt | xargs -n 1  ./ann_cons.py 
./comb_cons.py
./comp_cons.py
./merged_cons.py
./network_stat.py cons
./network_stat.py cons_classA
./statistics.py clean_consensus
./heatmaps_Gs_vs_Gio_cons.py
./heatmap.py cons

