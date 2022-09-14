1.network_script
contacts_cons.sh # Order of scripts to run
3dcontacts.sh #Calculate contacts between GPCR and Gprotein in pdb structures
contacts_cons.py #Map contacts to positions from sifts for each of pdb
ann_cons.py #Map contacts to GPCRDb positions
comb_cons.py #Combine all contact files to one
comp_cons.py #Calculate contact conservation grouped on gprotein family and calculate links and nodes for each gpcr gprotein network(all parameters for network)
consensus_struct.py #Calculate which pdb structures are representative based on smallest resolution and better coverage of Gprotein sequence
merged_cons.py #Creates final links and nodes file for cytoscape
network_stat.py #Statistics of network (node degree)
statistics.py #Statistics of pdb complexes(number,gprotein,family ,gpcr class)
centrality.py #Plots for centrality betweenes of networks

2. heatmap_script contains following scripts
barplot_cluster.py #Script that creates barplots of GPCR in different cluster and their gprotein coupling(UCM coupling)
heatmaps_Gs_vs_Gio_cons.py #Calculates log odds ratio for gpcr and gprotein contacts seperately and plots heatmaps of the contacts for all GPCRs from pdb complexes (all contacts,class A only and 0.2 cuttoff)
heatmap.py #Plot heatmaps with contact pairs for each gprotein coupling group with fraction of conservation(cuttoff 0.2)
heatmap_all_pair.py #Plots heatmaps of the contacts for all GPCRs from pdb complexes (0.1 cuttoff)
heatmap_all_pair_logodd.py #Calculates log odds ratio for gpcr and gprotein pair contact and plots heatmaps of the contacts for all GPCRs from pdb complexes (0.1 cuttoff conservation and >abs(2) log odds ratio)
heatmaps_Gs_vs_Gio_all1.py #Plots heatmaps of the contacts for all GPCRs from pdb complexes (without log odds)
logs.py # Plot distribution of log odds ratio for gpcr and gprotein 
logs_pair.py  #Plot distribution of log odds ratio for pair positions

3. rmsd_script contains scripts for calculation of rmsd
rmsd.sh #scripts to run to create rmsd files to calculate all rmsd between possible pairs of structures
Uniprot2PDBres_via_SIFTsxml.py #map uniprot and pdb positions via SIFT
bio_align_fit_n_rmsd.py #calculate rmsd and rmsf of pair of structures
extract_rmsd.py #combine all rmsd to one file
extract_rmsf.py  #combine all rmsf to one file
gen_bioalign.py #create a bash script with commands to run rmsd calculations for all pairs
rmsd_plot.py #plot rmsds and calculate statistics
rmsd_pos.py # calculate which gpcr positions and gprotein  consensus positions to use in to calculate rmsd (positions that show in each gpcr or gprotein of the structures)
rmsf_plot.py #plot rmsf for each of gprotein family for structure that is best fitting for the gprotein family
4. use_file #Annotation files
