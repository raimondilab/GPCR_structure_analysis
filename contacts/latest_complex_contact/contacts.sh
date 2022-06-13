#!/bin/bash
awk -F '\t' '{print $1}' ./use_file/GPCR_structs_clean.tsv> pdb_list.txt
for pdb in `cat pdb_list.txt`; do echo "/data/DB/pdb-mmCIF/"${pdb:1:2}"/"$pdb".cif.gz" >> pdb_list_new.txt ; done 
sed -i '1d' pdb_list_new.txt
xargs -a pdb_list_new.txt cp -t ./cifs/
./3dcontacts.sh 
ls ./cifs/*contacts.txt | xargs -n 1  ./contacts.py
ls ./cont_file/*cont.txt | xargs -n 1  ./ann.py 
./comb.py
./comp.py
./merged.py
./network_stat.py all
./network_stat.py all_classA
./statistics.py clean
./heatmaps_Gs_vs_Gio_all.py
./heatmap.py all
