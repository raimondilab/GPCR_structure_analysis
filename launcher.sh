#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=168:00:00
#PBS -q q07anacreon
#PBS -N GPCR_update

cd /home/pmiglionico/GPCR_structure_analysis
bash update_list.sh
bash update_stats.sh
