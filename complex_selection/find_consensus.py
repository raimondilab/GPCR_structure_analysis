#!/usr/bin/python3

'''This module generates 3 figures for the interface analysis. One using all the structures, one using all
representative structures and one averaging all the structures representing the same GPCR-Galpha pair'''

import pandas as pd

total = pd.read_csv("GPCR_structs_clean.tsv", sep='\t')

# Select representative structures: choose the one with the best resolution,
# in case of parity the one with the highest sequence coverage.
representative = total.sort_values(by=["Resolution", "Gprotein_residues", "Receptor_residues"], \
     key=lambda x: 10**9*total["Resolution"]-10**4*total["Receptor_residues"]-total["Gprotein_residues"])
representative = representative.drop_duplicates(subset=['Receptor_Gene_name', 'Gprotein_Gene_name'])
representative.to_csv("GPCR_consensus.tsv", sep='\t', index=False)
