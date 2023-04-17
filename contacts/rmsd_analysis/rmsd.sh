#!/bin/bash
#PBS -l select=1:ncpus=1 
#PBS -q q07anacreon
#PBS -N rmsd
source activate /home/mmatic/.conda/envs/transformer
cd  /home/mmatic/GPCR_structure_analysis/rmsd_all
##### First part of rmsd analysis exctracting position of GPCR and Gprotein of the structures to do rmsd calculation and creating files to run it
python Uniprot2PDBres_via_SIFTsxml.py
python rmsd_pos.py
python gen_bioalign.py
chmod 777 /home/mmatic/GPCR_structure_analysis/rmsd_all/rmsd_sh/*.sh
for file in /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all/*.sh; do qsub $file; done
