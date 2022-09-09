#!/data/SW/anaconda3/envs/myenv/bin/python

import sys
from filesplit.split import Split

fs="\t"
path='./cifs/'
Gprot_r={}
Gpcr_r={}
for l in open("Gprot_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split('_')[0]
  chain=l.split(fs)[0].split('_')[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("gpcr_pos_cons.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_cons_gpcr_cons.sh", "a") as f:
    print ("#!/bin/bash",file=f)
    for ii in range(len(Gpcr_r.keys())):
        for jj in range(ii+1, len(Gpcr_r.keys())):
            pdb1=list(Gpcr_r.keys())[ii]
            rc1=Gpcr_r[pdb1][0]
            r_res1=Gpcr_r[pdb1][1]
            pdb2=list(Gpcr_r.keys())[jj]
            rc2=Gpcr_r[pdb2][0]
            r_res2=Gpcr_r[pdb2][1]            
            mc1=Gprot_r[pdb1][0]
            m_res1=Gprot_r[pdb1][1].strip('"')
            mc2=Gprot_r[pdb2][0]
            m_res2=Gprot_r[pdb2][1].strip('"')
            print ("./bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif',path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f)
###
Gprot_r={}
Gpcr_r={}
for l in open("gprot_pos_cons.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("gpcr_pos_cons.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_all_gpcr_cons.sh", "a") as f1:
    print ("#!/bin/bash",file=f1)
    for ii in range(len(Gpcr_r.keys())):
        for jj in range(ii+1, len(Gpcr_r.keys())):
            pdb1=list(Gpcr_r.keys())[ii]
            rc1=Gpcr_r[pdb1][0]
            r_res1=Gpcr_r[pdb1][1]
            pdb2=list(Gpcr_r.keys())[jj]
            rc2=Gpcr_r[pdb2][0]
            r_res2=Gpcr_r[pdb2][1]            
            mc1=Gprot_r[pdb1][0]
            m_res1=Gprot_r[pdb1][1].strip('"')
            mc2=Gprot_r[pdb2][0]
            m_res2=Gprot_r[pdb2][1].strip('"')
            print ("./bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f1)

###
Gprot_r={}
Gpcr_r={}
for l in open("Gprot_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split('_')[0]
  chain=l.split(fs)[0].split('_')[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("gpcr_pos_cons_12.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_cons_gpcr_cons12.sh", "a") as f2:
    print ("#!/bin/bash",file=f2)
    for ii in range(len(Gpcr_r.keys())):
        for jj in range(ii+1, len(Gpcr_r.keys())):
            pdb1=list(Gpcr_r.keys())[ii]
            rc1=Gpcr_r[pdb1][0]
            r_res1=Gpcr_r[pdb1][1]
            pdb2=list(Gpcr_r.keys())[jj]
            rc2=Gpcr_r[pdb2][0]
            r_res2=Gpcr_r[pdb2][1]            
            mc1=Gprot_r[pdb1][0]
            m_res1=Gprot_r[pdb1][1].strip('"')
            mc2=Gprot_r[pdb2][0]
            m_res2=Gprot_r[pdb2][1].strip('"')
            print ("./bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f2)

###
Gprot_r={}
Gpcr_r={}
for l in open("gprot_pos_cons.tsv", "rt"):
    pdb=l.split(fs)[0]
    chain=l.split(fs)[1]
    res=l.split(fs)[2].strip("\n").strip(",")
    #print (pdb, len(res.split(",")))
    Gprot_r[pdb]=[chain,res]
for l in open("gpcr_pos_cons_12.tsv", "rt"):
    pdb=l.split(fs)[0]
    chain=l.split(fs)[1]
    res=l.split(fs)[2].strip("\n").strip(",")
    #print (pdb, len(res.split(",")))
    Gpcr_r[pdb]=[chain,res]

#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_all_gpcr_cons12.sh", "a") as f3:
    print ("#!/bin/bash",file=f3)
    for ii in range(len(Gpcr_r.keys())):
        for jj in range(ii+1, len(Gpcr_r.keys())):
            pdb1=list(Gpcr_r.keys())[ii]
            rc1=Gpcr_r[pdb1][0]
            r_res1=Gpcr_r[pdb1][1]
            pdb2=list(Gpcr_r.keys())[jj]
            rc2=Gpcr_r[pdb2][0]
            r_res2=Gpcr_r[pdb2][1]            
            mc1=Gprot_r[pdb1][0]
            m_res1=Gprot_r[pdb1][1].strip('"')
            mc2=Gprot_r[pdb2][0]
            m_res2=Gprot_r[pdb2][1].strip('"')
            print ("./bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f3)

split0 = Split('./rmsd_sh/bioalign_gprot_cons_gpcr_cons12.sh','./rmsd_sh/')
split1 = Split('./rmsd_sh/bioalign_gprot_cons_gpcr_cons.sh','./rmsd_sh/')
split2 = Split('./rmsd_sh/bioalign_gprot_all_gpcr_cons.sh','./rmsd_sh/')
split3 = Split('./rmsd_sh/bioalign_gprot_all_gpcr_cons12.sh','./rmsd_sh/')

split0.bylinecount(linecount=2000, includeheader = True)
split1.bylinecount(linecount=2000, includeheader = True) 
split2.bylinecount(linecount=2000, includeheader = True) 
split3.bylinecount(linecount=2000, includeheader = True)  
