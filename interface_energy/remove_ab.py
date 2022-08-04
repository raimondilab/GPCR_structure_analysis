#!/usr/bin/env python3

'''This module finds the chains corresponding to antibodies or nanobodies
in a .cif structure and removes them saving the structure in .pdb format'''

import sys
import warnings
import Bio.PDB

warnings.filterwarnings("ignore")
Ab = []

class ResSelect(Bio.PDB.Select):
    def accept_residue(self, res):
        return res.parent.id not in Ab

pdb = sys.argv[1].split('/')[-1][:4]
print(pdb)

structdict = Bio.PDB.MMCIF2Dict.MMCIF2Dict(sys.argv[1])
for i in range(len(structdict['_entity.pdbx_description'])):
    descr = str.lower(structdict['_entity.pdbx_description'][i])
    # This list of if statements checks if in the entity desctiption there is the word "Antibody",
    # "Nanobody", "Fab" or "scF", it also checks for some typos that can be found here and there
    if descr.find('body') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)
    elif descr.find('scf') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)
    elif descr.find('fab') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)
    elif descr.find('svf') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)
    elif descr.find('nb') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)
    elif descr.find('boy') != -1:
        print(structdict['_entity.pdbx_description'][i])
        Ab.append(i)

# Useful to check if the antibody is missing or if there are more words to add to the previous loop
if Ab == []:
    print(structdict['_entity.pdbx_description'])

# Convert the entity number to chain_id
for i in range(len(Ab)):
    Ab[i] = structdict['_entity_poly.pdbx_strand_id'][Ab[i]]

# Upload the structure and save it without antibodies
struct = Bio.PDB.MMCIFParser().get_structure(pdb, sys.argv[1])
io = Bio.PDB.PDBIO()
io.set_structure(struct)
io.save(sys.argv[2], ResSelect())
