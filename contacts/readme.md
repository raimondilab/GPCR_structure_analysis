2. heatmap_script contains following scripts
barplot_cluster.py
#Script that creates barplots of GPCR in different cluster and their gprotein coupling(UCM coupling)
heatmaps_Gs_vs_Gio_cons.py
#Calculates log odds ratio for gpcr and gprotein contacts seperately and plots heatmaps of the contacts for all GPCRs from pdb complexes (all contacts,class A only and 0.2 cuttoff)
heatmap.py
#Plot heatmaps with contact pairs for each gprotein coupling group with fraction of conservation(cuttoff 0.2)
heatmap_all_pair.py
#Plots heatmaps of the contacts for all GPCRs from pdb complexes (0.1 cuttoff)
heatmap_all_pair_logodd.py
#Calculates log odds ratio for gpcr and gprotein pair contact and plots heatmaps of the contacts for all GPCRs from pdb complexes (0.1 cuttoff conservation and >abs(2) log odds ratio)
heatmaps_Gs_vs_Gio_all1.py
#Plots heatmaps of the contacts for all GPCRs from pdb complexes (without log odds)
logs.py
# Plot distribution of log odds ratio for gpcr and gprotein 
logs_pair.py
# Plot distribution of log odds ratio for pair positions
