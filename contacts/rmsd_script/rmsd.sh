#!/bin/bash
awk -F '\t' '{print $1}' GPCR_structs_clean_consensus.tsv > pdb_list.txt
sed -i '1d' pdb_list.txt
for pdb in `cat pdb_list.txt`; do echo "/data/DB/pdb-mmCIF/"${pdb:1:2}"/"$pdb".cif.gz" >> pdb_list_new.txt ; done
xargs -a pdb_list_new.txt cp -t ./cifs
printf '1\ni\n#PDBID\tCHAIN\tPDB_RES_NUM\tUNIPROT_AC\tUNIPROT_RES_NUM\n.\nw\n' | ed -s uniprot_pdb.tsv
xargs -a pdb_list.txt |xargs -n 1 ./Uniprot2PDBres_via_SIFTsxml.py >> uniprot_pdb.tsv
./gprot.sh
./rmsd_pos.py
sed -i '1d' gprot_pos_cons.tsv
sed -i '1d' gpcr_pos_cons.tsv
sed -i '1d' gpcr_pos_cons_12.tsv
./gen_bioalign.py
