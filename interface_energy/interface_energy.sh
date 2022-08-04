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
	if [ $pdb == "PDB_ID" ]
	then
		continue
	fi
	if [ -f relaxed_struct_no_Ab/${pdb}.pdb ]
	then
		continue
	fi
    if ((n < 52))
    then
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/relax.default.linuxgccrelease -s structures_no_Ab/$pdb\.pdb -ignore_unrecognized_res -remember_unrecognized_res -out:path:pdb relaxed_struct_no_Ab -out:file:scorefile relaxed_struct_no_Ab/$pdb\_scorefile.sc -nstruct 1 &
		n=$((n+1)) 
    else
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/relax.default.linuxgccrelease -s structures_no_Ab/$pdb\.pdb -ignore_unrecognized_res -remember_unrecognized_res -out:path:pdb relaxed_struct_no_Ab -out:file:scorefile relaxed_struct_no_Ab/$pdb\_scorefile.sc -nstruct 1
		n=0   
    fi
done < ../GPCR_structure_analysis/GPCR_structs_clean.tsv
wait

while read model
do
	pdb=`echo $model | awk -F ' ' '{print $1}'`
	if [ $pdb == "PDB_ID" ]
	then
		continue
	fi
	c1=`echo $model | awk -F ' ' '{print $5}'`
	c2=`echo $model | awk -F ' ' '{print $9}'`
	rm relaxed_struct_no_Ab/${pdb}_interface.tsv
	mv relaxed_struct_no_Ab/$pdb\_0001.pdb relaxed_struct_no_Ab/$pdb\.pdb
    if ((n < 52))
    then
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.linuxgccrelease -s relaxed_struct_no_Ab/$pdb\.pdb -interface $c1\_$c2 -out:file:score_only relaxed_struct_no_Ab/$pdb\_interface.tsv -pack_input -pack_separated -ignore_unrecognized_res &
		n=$((n+1)) 
    else
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/InterfaceAnalyzer.linuxgccrelease -s relaxed_struct_no_Ab/$pdb\.pdb -interface $c1\_$c2 -out:file:score_only relaxed_struct_no_Ab/$pdb\_interface.tsv -pack_input -pack_separated -ignore_unrecognized_res
		n=0   
    fi
done < ../GPCR_structure_analysis/GPCR_structs_clean.tsv
wait
