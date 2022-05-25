#!/bin/bash

mkdir ../GPCR_experimental_structures
qsub relax.sh
python3 find_interface_energy.py GPCR_structs.tsv ../GPCR_experimental_structures/binding_energy.tsv
mkdir ../GPCR_experimental_structures/figures
python3 interface_energy_plot.py
