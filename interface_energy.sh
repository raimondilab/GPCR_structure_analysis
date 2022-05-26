#!/bin/bash
#PBS -l select=1:ncpus=52
#PBS -l walltime=168:00:00
#PBS -q q07anacreon
#PBS -N interface_energy

module load gcc
n=0
cd "/home/pmiglionico/GPCR_experimental_structures"

while read model
do
	pdb=`echo $model | awk -F ' ' '{print $1}'`
	c1=`echo $model | awk -F ' ' '{print $5}'`
	c2=`echo $model | awk -F ' ' '{print $9}'`
	if [ -f "relaxed_struct/${pdb}_interface.tsv" ]
	then
		continue
	fi
	mv relaxed_struct/$pdb\_0001.pdb relaxed_struct/$pdb\.pdb
    if ((n < 52))
    then
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.linuxgccrelease -s relaxed_struct/$pdb\.pdb -interface $c1\_$c2 -out:file:score_only relaxed_struct/$pdb\_interface.tsv -pack_input -pack_separated -ignore_unrecognized_res &
		n=$((n+1)) 
    else
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.linuxgccrelease -s relaxed_struct/$pdb\.pdb -interface $c1\_$c2 -out:file:score_only relaxed_struct/$pdb\_interface.tsv -pack_input -pack_separated -ignore_unrecognized_res
		n=0   
    fi
done < ../GPCR_analysis/GPCR_structs_clean.tsv
wait
