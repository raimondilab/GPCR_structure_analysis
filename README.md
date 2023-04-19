# GPCR_structure_analysis

This repository contains the code needed to repeat the structural analysis of the GPCR-Galpha experimental structures. The paths are set to make the scripts run on our cluster, if you would like to run it on your machine you need to change the paths of the databases. Feel free to reach out for any question.

To update the list of structures and run the analysis, run:
```
bash update_list.sh
```

To repeat the analysis on the new list, run:
```
bash update_stats.sh
```

This command is quite costly (about 6 hours of CPU time for each new structure), because it needs to relax the structures before computing the interface energy.

# pymol_chimera_sessions
Folder containg for each figure a session either in pymol or chimerax
# contacts
Folder containg contact analysis and rmsd analysis scripts
