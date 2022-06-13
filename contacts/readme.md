1. latest_complex_contact contain all the complex structure  and their contacts and plots refresh GPCR_struct_clean.tsv and run contacts.sh
 2. consensus contain consensus   complex structure  and their contacts and plots refresh GPCR_struct_clean.tsv and run contacts_cons.sh
3. rmsd contains rmsd calculation and all the plots  refresh GPCR_struct_clean.tsv
4. run /rmsd.sh  then ex_rmsd.sh (it takes a while it calculates rmsd (it sets for parallel running) for 4 different combinations of gpcr and gprot) and then plot.sh to get plots 
5. rmsd_consensus contains rmsd calculation and all the plots for consensus sequences refresh GPCR_struct_clean.ts)
6. run /rmsd.sh  then ex_rmsd.sh (it takes a while it calculates rmsd for 4 different combinations of gpcr and gprot) and then plot.sh to get plots it recalculates everything for the consensus set
 folders with rmsd calculation
 bioalign_gprot_all_gpcr_all (all the common positions for gpcr and gprotein for all the structures positons based on pdb to uniprot mapping) 
 bioalign_gprot_all_gpcr_all12 (all the common positions for gprot and 12  common in all for gpcr(1.49,1.5,1.51,3.49,3.5,3.51,4.49,4.5,4.51,5.41,6.49,6.5) for all the structures based on pdb to uniprot mapping) 
 bioalign_gprot_cons_gpcr_all(all the common positions for gprot based on msa alignment and  all in gpcr based on pdb to uniprot mapping) 
 bioalign_gprot_cons_gpcr_all(all the common positions for gprot based on msa alignment and 12  common in all for gpcr(1.49,1.5,1.51,3.49,3.5,3.51,4.49,4.5,4.51,5.41,6.49,6.5) for all the structures based on pdb to uniprot mapping) 
 same in rmsd_consensus except consensus structures and added one more position 6.51
 
 excluded 5g53 cause it fails running the selection part
 excluded 7jjo,6p9x,6pb0 after rmsd (in plot) cause of too big values
