#!/data/SW/anaconda3/envs/myenv/bin/python

"""
Sequence-based structural alignment of two proteins.
"""

from __future__ import print_function, division

import argparse
import os, gzip
import numpy as np
from Bio.PDB import PDBParser, FastMMCIFParser, Superimposer, PDBIO, QCPSuperimposer
from Bio.PDB.Polypeptide import is_aa

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1



def rmsd(sele1, sele2):
###Adapted from SVDSuperimposer module in Bio.PDB
    """Return rms deviations between coords1 and coords2 (PRIVATE)."""
    coords1=[atom.coord for atom in sele1]
    coords2=[atom.coord for atom in sele2]
    coords1=np.array(coords1)
    coords2=np.array(coords2)    
    diff = coords1 - coords2
    return np.sqrt(sum(sum(diff * diff)) / coords1.shape[0])


def parse_structure(spath):
    """Parses a PDB/cif structure"""

    if not os.path.isfile(spath):
        return IOError('File not found: {0}'.format(spath))

    if spath.endswith(('pdb', 'ent')):
        parser = PDBParser()
    elif spath.endswith('cif'):
        parser = FastMMCIFParser()
    else:
        raise Exception('Format not supported ({0}). Must be .pdb/.ent or .cif'.format(spath))

    sname = os.path.basename(spath.split('.')[0])
    return parser.get_structure(sname, spath)


if __name__ == '__main__':

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-refe', help='Reference Structure')
    ap.add_argument('-mobi', help='Mobile Structure')
    ap.add_argument('-r_fit_chain', default='A', help='Reference Structure Chain for fitting')
    ap.add_argument('-m_fit_chain', default='A', help='Mobile Structure Chain for fitting')
    ap.add_argument('-r_fit_sele',  default='A', help='Reference Structure Selection for fitting')
    ap.add_argument('-m_fit_sele',  default='A', help='Mobile Structure Selection for fitting')
    ap.add_argument('-r_rms_chain', default='A', help='Reference Structure Chain for RMSD')
    ap.add_argument('-m_rms_chain', default='A', help='Mobile Structure Chain for RMSD')
    ap.add_argument('-r_rms_sele',  default='A', help='Reference Structure Selection for RMSD')
    ap.add_argument('-m_rms_sele',  default='A', help='Mobile Structure Selection for RMSD')

    cline = ap.parse_args()

    # Parse structures; we keep all the chain as we frequently want to calculate the RMSD on a selection different from the one used to perform the fitting
    reference = parse_structure(cline.refe)
    mobile = parse_structure(cline.mobi)

    print (len(cline.r_fit_sele.split(",")), len(cline.m_fit_sele.split(",")))

#####Atom lists for fitting
    refe_ca_list, mobi_ca_list = [], []
    for ii in range(len(cline.r_fit_sele.split(","))):
        refe_res=int(cline.r_fit_sele.split(",")[ii])
        mobi_res=int(cline.m_fit_sele.split(",")[ii])
        refe_ca_list.append(reference[0][cline.r_fit_chain][refe_res]['CA'])
        mobi_ca_list.append(mobile[0][cline.m_fit_chain][mobi_res]['CA'])


    # Superimpose matching residues
    si = Superimposer()
    si.set_atoms(refe_ca_list, mobi_ca_list)

    # Transform & Write Mobile
    si.apply(mobile.get_atoms())
    print ("Fitting RMSD = ", si.rms)

###After superposition using the 1st selection, we calculate the RMSD on the second selection
####Atom lists for RMSD in place (we don't fit another time)
    refe_ca_list2, mobi_ca_list2 = [], []
    for ii in range(len(cline.r_rms_sele.split(","))):
        refe_res2=int(cline.r_rms_sele.split(",")[ii])
        mobi_res2=int(cline.m_rms_sele.split(",")[ii])
        refe_ca_list2.append(reference[0][cline.r_rms_chain][refe_res2]['CA'])
        mobi_ca_list2.append(mobile[0][cline.m_rms_chain][mobi_res2]['CA'])
    
    #si2 = Superimposer()
    #si2.set_atoms(refe_ca_list2, mobi_ca_list2)
    #print ("Selection RMSD = ", si2.rms)
    print ("Selection RMSD = ", rmsd(refe_ca_list2, mobi_ca_list2)) 

    io = PDBIO()
    io.set_structure(mobile)
    m_transformed_name = '%s_%s_transformed.pdb' % ((cline.refe.split(".cif")[0]).split("/")[-1], (cline.mobi.split(".cif")[0]).split("/")[-1])
    io.save(m_transformed_name)
