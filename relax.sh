#!/bin/bash
#PBS -l select=1:ncpus=52
#PBS -l walltime=168:00:00
#PBS -q q07anacreon
#PBS -N relax

module load gcc
n=0
cd "/home/pmiglionico/GPCR_experimental_structures"
mkdir structures
mkdir relaxed_struct

while read model
do
	pdb=`echo $model | awk -F ' ' '{print $1}'`
	if [ ! -f "structures/${pdb}.cif" ]
	then
		cp /home/fraimondi/BIOINFO1_DB/pdb-mmCIF/${pdb:1:2}/$pdb\.cif.gz structures/
		gunzip structures/$pdb\.cif.gz
	fi
	if [ ! -f "/home/pmiglionico/pdb-mmCIF_CBcontacts_blast/${pdb:1:2}/${pdb}_3dc.txt" ]
	then
		/home/pmiglionico/scripts/cifparse-obj-v7.105-prod-src/parser-test-app/bin/3DContact_nodom structures/${pdb}.cif 8 > /home/pmiglionico/pdb-mmCIF_CBcontacts_blast/${pdb:1:2}/${pdb}_3dc.txt &
	fi
	if [ -f "relaxed_struct/${pdb}.pdb" ]
	then
		continue
	fi
    if ((n < 52))
    then
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/relax.default.linuxgccrelease -s structures/$pdb\.cif -ignore_unrecognized_res -remember_unrecognized_res -out:path:pdb relaxed_struct -out:file:scorefile relaxed_struct/$pdb\_scorefile.sc -nstruct 1 &
		n=$((n+1)) 
    else
		/home/pmiglionico/rosetta_src_2021.16.61629_bundle/main/source/bin/relax.default.linuxgccrelease -s structures/$pdb\.cif -ignore_unrecognized_res -remember_unrecognized_res -out:path:pdb relaxed_struct -out:file:scorefile relaxed_struct/$pdb\_scorefile.sc -nstruct 1
		n=0   
    fi
done < ../GPCR_analysis/GPCR_structs.tsv
wait
