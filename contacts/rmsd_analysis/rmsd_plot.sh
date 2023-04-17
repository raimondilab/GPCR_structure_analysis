#!/bin/bash
#PBS -l select=1:ncpus=1 
#PBS -q q07anacreon
#PBS -N rmsd
source activate /home/mmatic/.conda/envs/transformer
cd  /home/mmatic/GPCR_structure_analysis/rmsd_all
#####Extract the rmsd and rmsf after calculations have finished and plot rmsd and rmsf

python extract_rmsd.py 
python extract_rmsd.py /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all/ > ./plots/bioalign_gprot_all_gpcr_all.txt
python rmsd_plot.py bioalign_gprot_all_gpcr_all
python extract_rmsf.py bioalign_gprot_all_gpcr_all
python rmsf_plot.py gprot_all_gpcr_all
