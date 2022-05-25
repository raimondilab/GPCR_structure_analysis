# GPCR_structure_analysis

This repository contains the code needed to repeat the structural analysis of the GPCR-Galpha experimental structures. The paths are set to make the scripts run on the trantor cluster.

To update the list of structures to analyze, run:
```
bash update_list.sh
```
The output is a file named GPCR_structs.tsv. This file needs to be manually reviewed, in particular to remove duplicated Galpha-GPCR pairs that are due to dimeric or chimeric complexes.

After the revision, to update all the analysis with the new set of structures, run:
```
bash run_update.sh
```
This command is quite costly (about 6 hours of CPU time for each new structure), because it needs to relax the structures before computing the interface energy.
