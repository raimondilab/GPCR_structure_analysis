#!/bin/bash

source activate bioinfo

# Find a list of pdb structures
echo "Generating list of structures..."
zgrep PF00503 /home/fraimondi/BIOINFO1_DB/SIFTS/pdb_chain_pfam.tsv.gz > Galpha_list.txt
python3 find_GPCR.py
rm Galpha_list.txt

# Compute relaxed structures and contacts (if not already available)
echo "Relaxing structures..."
mkdir ../GPCR_experimental_structures
qsub -sync y relax.sh

# Resolve ambiguities in GPCR or Gproteins identifiers or chains in the structure and add informations to the sturcture file
echo "Cleaning list of pairs..."
python3 clean_selection.py

# Compute interface energy
echo "Computing Interface energies..."
qsub -sync y interface_energy.sh
python3 find_interface_energy.py GPCR_structs_clean.tsv ../GPCR_experimental_structures/binding_energy.tsv
echo "Computed interface energy"

# Draw figures
echo "Drawing figures..."
mkdir ../GPCR_experimental_structures/figures
python3 interface_energy_plot.py
