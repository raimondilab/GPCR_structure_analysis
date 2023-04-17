# Contact analysis 
This folder contains scripts to run contact analysis from the paper.
### contacts.sh
 Running contact analysis order
### copy_cifs.py
  Copy cifs files, extracted contact interfaces for each cif file into personal folder
### ann.py
  Annotate contact interface for each structure with GPCRDB numbering and CGN numbering
### comb.py
   Combine all interface files into one removing redundant interface contacts
### comp.py
   Calculate contact fractions for each G protein family and network link(edge) and node files
### merged.py
   Combine G protein family links and nodes files into master links and node files to use for network   building in Cytoscape
### network_stat.py
   Calculate network statistics degree, radius, and centrality
### statistics.py
   Create a sunburst diagram of the experimental set of structures annotating with GPCR class
### logs.py
   Calculate and plot log odds ratio of Gs vs Gio structure contacts separately for GPCR and G protein portion
### logs_pair.py
   Calculate and plot log odds ratio of Gs vs Gio structure contacts for paired GPCR-G protein contacts
### heatmap.py
   Plot heatmaps of Gs vs Gio contacts shown in at least 20% of GPCRs and conservation of those contacts annotated with log odds ratio
### heatmaps_Gs_vs_Gio_all1.py
   Plot heatmaps with different filters of contact pairs or GPCR and G protein separately for all structures
### heatmaps_Gs_vs_Gio_logodd.py
   Plot heatmaps with different filters of contact pairs or GPCR and G protein separately for Gs vs Gio structures and annotated with log odds ratio

Folder ann_file contains each of the contact interface for experimental structures annotated with GPCRDB numbering and CGN numbering. 

