# GPCR_structure_analysis

This repository contains the code needed to repeat the structural analysis of the GPCR-Galpha experimental structures. The paths are set to make the scripts run on the trantor cluster.

To update the list of structures and run the analysis, run:
```
bash run_update.sh
```

This command is quite costly (about 6 hours of CPU time for each new structure), because it needs to relax the structures before computing the interface energy.
