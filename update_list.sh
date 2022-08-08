#!/bin/bash

source activate bioinfo

# Find a list of pdb structures
echo "Generating list of structures..."
zgrep PF00503 /home/pmiglionico/pdb_chain_pfam.tsv.gz > Galpha_list.txt
python3 complex_selection/find_GPCR.py
rm Galpha_list.txt

# Copy structures to working folder
mkdir -p ../GPCR_experimental_structures/structures
while read model
do
	pdb=`echo $model | awk -F ' ' '{print $1}'`
	if [ $pdb == "PDB_ID" ]
	then
		continue
	fi
	if [ ! -f "../GPCR_experimental_structures/structures/${pdb}.cif" ]
	then
		cp /home/fraimondi/BIOINFO1_DB/pdb-mmCIF/${pdb:1:2}/$pdb\.cif.gz ../GPCR_experimental_structures/structures/
		gunzip ../GPCR_experimental_structures/structures/$pdb\.cif.gz
	fi
done < GPCR_structs.tsv

# Resolve ambiguities in GPCR or Gproteins identifiers or chains in the structure and add informations to the structure file
echo "Cleaning list of pairs..."
python3 complex_selection/clean_selection.py
mv GPCR_structs_clean.tsv GPCR_structs.tsv
mkdir ../GPCR_experimental_structures/cont_file
python3 complex_selection/check_contacts.py
rm GPCR_structs.tsv
python3 complex_selection/find_consensus.py
