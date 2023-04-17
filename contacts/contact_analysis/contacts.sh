#!/bin/bash
#PBS -l select=1:ncpus=5 
#PBS -q q07anacreon
#PBS -N STRUCTURE
source activate /home/mmatic/.conda/envs/transformer
cd /home/mmatic/GPCR_structure_analysis/contact_all
#################### Running contact analysis order
python copy_cifs.py
python ann.py
python comb.py
python comp.py
python merged.py
python network_stat.py all 
python network_stat.py all_classA
python statistics.py clean
python logs.py
python logs_pair.py
python heatmaps_Gs_vs_Gio_all1.py
python heatmap.py all
python heatmaps_Gs_vs_Gio_logodd.py
