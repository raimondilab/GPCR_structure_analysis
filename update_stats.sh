#!/bin/bash

source activate bioinfo

# Remove the antibody chain from the structures
for pdb in `ls -1t ../GPCR_experimental_structures/structures/*cif`
do
    pdbout=`echo $pdb | cut -f 4 -d "/"`
    if [ -f ../GPCR_experimental_structures/structures_no_Ab/${pdbout%%.cif}.pdb ]
	then
		continue
	fi
	python3 interface_energy/remove_ab.py $pdb ../GPCR_experimental_structures/structures_no_Ab/${pdbout%%.cif}.pdb
done

echo "Relaxing structures and computing interface energies..."
#qsub -Wblock=true interface_energy/interface_energy.sh

# Add the interface energies to the table
python3 interface_energy/find_interface_energy.py GPCR_structs_clean.tsv binding_energy.tsv
mv binding_energy.tsv GPCR_structs_clean.tsv

# Find the best structure for each GPCR-Gprotein complex
python3 complex_selection/find_consensus.py
python3 complex_selection/isrepresentative.py GPCR_consensus.tsv GPCR_structs_clean.tsv GPCR_isrepresentative.tsv
mv GPCR_isrepresentative.tsv GPCR_structs_clean.tsv

# Draw figures
echo "Drawing figures..."
mkdir ../GPCR_experimental_structures/figures
python3 interface_energy/interface_energy_plot.py
python3 interface_energy/interface_energy_plot_classA.py
