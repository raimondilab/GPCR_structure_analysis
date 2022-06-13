#!/bin/bash
awk -F '\t' '{print $1}' GPCR_structs_clean.tsv> pdb_list.txt
sed -i '1d' pdb_list.txt
for pdb in `cat pdb_list.txt`; do echo "/data/DB/pdb-mmCIF/"${pdb:1:2}"/"$pdb".cif.gz" >> pdb_list_new.txt ; done
xargs -a pdb_list_new.txt cp -t ./cifs
printf '1\ni\n#PDBID\tCHAIN\tPDB_RES_NUM\tUNIPROT_AC\tUNIPROT_RES_NUM\n.\nw\n' | ed -s uniprot_pdb.tsv
xargs -a pdb_list.txt |xargs -n 1 ./Uniprot2PDBres_via_SIFTsxml.py >> uniprot_pdb.tsv
./gprot.sh
./rmsd_pos.py
sed -i '1d' gprot_pos_all.tsv
sed -i '1d' gpcr_pos_all.tsv
sed -i '1d' gpcr_pos_all_12.tsv
./gen_bioalign.py
chmod 777 ./rmsd_sh/*.sh
for((i=1;i<17;i++)); do nohup bash ./rmsd_sh/bioalign_gprot_all_gpcr_all_${i}.sh & done
mv *fit.txt ./bioalign_gprot_all_gpcr_all/
mv *.pdb ./bioalign_gprot_all_gpcr_all
for((i=1;i<17;i++)); do nohup bash ./rmsd_sh/bioalign_gprot_cons_gpcr_all_${i}.sh & done
mv *fit.txt ./bioalign_gprot_cons_gpcr_all/
mv *.pdb ./bioalign_gprot_cons_gpcr_all/
for((i=1;i<17;i++)); do nohup bash ./rmsd_sh/bioalign_gprot_cons_gpcr_all12_${i}.sh & done
mv *fit.txt ./bioalign_gprot_cons_gpcr_all12/
mv *.pdb ./bioalign_gprot_cons_gpcr_all12/
for((i=1;i<17;i++)); do nohup bash ./rmsd_sh/bioalign_gprot_cons_gpcr_all12_${i}.sh & done
mv *fit.txt ./bioalign_gprot_all_gpcr_all12/
mv *.pdb ./bioalign_gprot_all_gpcr_all12/
./extract_rmsd.py	./bioalign_gprot_all_gpcr_all > rmsd_gprot_all_gpcr_all.txt
./extract_rmsd.py	./bioalign_gprot_cons_gpcr_all/ > rmsd_gprot_cons_gpcr_all.txt
./extract_rmsd.py ./bioalign_gprot_cons_gpcr_all12/ > rmsd_gprot_cons_gpcr_all12.txt
./extract_rmsd.py ./bioalign_gprot_all_gpcr_all12/ > rmsd_gprot_all_gpcr_all12.txt
