#!/bin/bash

source activate bioinfo

for pdb in `ls -1t ../GPCR_experimental_structures/structures/*cif`
do
    pdbout=`echo $pdb | cut -f 4 -d "/"`
    if [ -f ../GPCR_experimental_structures/structures_no_Ab/${pdbout%%.cif}.pdb ]
	then
		continue
	fi
	python3 interface_energy/remove_ab.py $pdb ../GPCR_experimental_structures/structures_no_Ab/${pdbout%%.cif}.pdb
done

# Compute relaxed structures and contacts (if not already available)
echo "Relaxing structures and computing interface energies..."
qsub -Wblock=true interface_energy/interface_energy.sh
python3 interface_energy/find_interface_energy.py GPCR_structs_clean.tsv ../GPCR_experimental_structures/binding_energy.tsv

# Draw figures
echo "Drawing figures..."
mkdir ../GPCR_experimental_structures/figures
python3 interface_energy/interface_energy_plot.py > ../GPCR_experimental_structures/stats.txt
python3 interface_energy/interface_energy_plot_classA.py > ../GPCR_experimental_structures/stats_classA.txt
