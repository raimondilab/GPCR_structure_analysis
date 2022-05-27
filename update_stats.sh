#!/bin/bash

source activate bioinfo

# Compute relaxed structures and contacts (if not already available)
echo "Relaxing structures and computing interface energies..."
#qsub -Wblock=true interface_energy.sh
python3 find_interface_energy.py GPCR_structs_clean.tsv ../GPCR_experimental_structures/binding_energy.tsv

# Draw figures
echo "Drawing figures..."
mkdir ../GPCR_experimental_structures/figures
python3 interface_energy_plot.py > ../GPCR_experimental_structures/stats.txt
