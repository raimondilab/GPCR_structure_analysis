# RMSD analysis
### rmsd.sh
First part of rmsd analysis exctracting position of GPCR and Gprotein of the structures to do rmsd(root mean square  deviation) and rmsf (root mean square fluctuation) calculation and creating files to run it
### Uniprot2PDBres_via_SIFTsxml.py
Mapping pdb positions in structure to uniprot positions in sequence for experimental structures
### rmsd_pos.py
Extracting common positions for rmsd calculation for GPCR and G protein part of the structures
### gen_bioalign.py
Creating bash files to run parallely multiple rmsd calculations for each of the pair of structures
### bio_align_fit_n_rmsd.py
Sequence-based structural alignment of two proteins and RMSD calculation
### rmsd_plot.sh
Extract the rmsd and rmsf after calculations have finished and plot rmsd and rmsf
### extract_rmsd.py
Extract rmsd from the calculated files
### extract_rmsf.py
Extract rmsf from the files
### rmsd_plot.py
Plot Heatmaps of rmsd calculation and the violin plots of the distribution of the rmsd.
### rmsf_plot.py
Plotting RMSF distribution of Gs vs Gio
