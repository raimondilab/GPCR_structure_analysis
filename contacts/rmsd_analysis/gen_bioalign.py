#!/data/SW/anaconda3/envs/myenv/bin/python

import sys
import os
######### Creating bash files to run parallely multiple rmsd calculations for each of the pair of structures

fs="\t"
path='/home/mmatic/GPCR_structure_analysis/rmsd_all/cifs/'
header='#!/bin/bash\n#PBS -l select=1:ncpus=1 \n#PBS -q q07anacreon\n#PBS -N rmsd\nsource activate /home/mmatic/.conda/envs/transformer\n'
#header='#!/bin/bash'
Gprot_r={}
Gpcr_r={}
for l in open("Gprot_pdb_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split('_')[0]
  chain=l.split(fs)[0].split('_')[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("./pos/gpcr_pos_all.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_cons_gpcr_all.sh", "w") as f:
    print (header + 'cd /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_cons_gpcr_all/\n',file=f)
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
            print ("python /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_cons_gpcr_all/bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif',path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f)
###
Gprot_r={}
Gpcr_r={}
for l in open("./pos/gprot_pos_all.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("./pos/gpcr_pos_all.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_all_gpcr_all.sh", "w") as f1:
    print (header + 'cd /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all/\n',file=f1)
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
            print ("python /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all/bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f1)

###
Gprot_r={}
Gpcr_r={}
for l in open("Gprot_pdb_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split('_')[0]
  chain=l.split(fs)[0].split('_')[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("./pos/gpcr_pos_all_12.tsv", "rt"):
  pdb=l.split(fs)[0]
  chain=l.split(fs)[1]
  res=l.split(fs)[2].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_cons_gpcr_all12.sh", "w") as f2:
    print (header+ 'cd /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_cons_gpcr_all12/\n',file=f2)
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
            print ("python /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_cons_gpcr_all12/bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f2)

###
Gprot_r={}
Gpcr_r={}
for l in open("./pos/gprot_pos_all.tsv", "rt"):
    pdb=l.split(fs)[0]
    chain=l.split(fs)[1]
    res=l.split(fs)[2].strip("\n").strip(",")
    #print (pdb, len(res.split(",")))
    Gprot_r[pdb]=[chain,res]
for l in open("./pos/gpcr_pos_all_12.tsv", "rt"):
    pdb=l.split(fs)[0]
    chain=l.split(fs)[1]
    res=l.split(fs)[2].strip("\n").strip(",")
    #print (pdb, len(res.split(",")))
    Gpcr_r[pdb]=[chain,res]

#print (list(Gpcr_r.keys()))
with open("./rmsd_sh/bioalign_gprot_all_gpcr_all12.sh", "w") as f3:
    print (header+ 'cd /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all12/\n',file=f3)
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
            print ("python /home/mmatic/GPCR_structure_analysis/rmsd_all/bioalign_gprot_all_gpcr_all12/bio_align_fit_n_rmsd.py -refe %s -mobi %s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (path+pdb1+'.cif', path+pdb2+'.cif', rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1, pdb2),file=f3)


def split_file(input_file,input_dir,output_dir,header_lines,num_files):
    input_file = input_file
    input_dir = input_dir
    output_dir=output_dir
    # Open the input file and read the lines
    with open(os.path.join(input_dir, input_file), "r") as f:
        lines = f.readlines()
        # Define the number of lines to keep in each output file
    header_lines = header_lines

    # Calculate the number of output files
    num_files = num_files

    # Calculate the number of lines to split into each output file
    num_lines = (len(lines) - header_lines) // num_files

    # Loop through the output files and write the lines
    for i in range(num_files):
        # Define the output file name
        output_file = f"{input_file[:-3]}_{i+1}.sh"
    
        # Open the output file for writing
        with open(os.path.join(output_dir, output_file), "w") as f:
            # Write the first 5 lines from the input file
            header_lines_list = lines[:header_lines]
            last_line =header_lines_list[header_lines-1].strip()+str(i+1)+'\n'
            f.writelines(header_lines_list[:header_lines-1])
            f.writelines(last_line)
        
            # Write the remaining lines, splitting them into equal-sized chunks
            start_idx = header_lines + i*num_lines
            end_idx = header_lines + (i+1)*num_lines if i < num_files-1 else None
            f.writelines(lines[start_idx:end_idx])


split_file('bioalign_gprot_cons_gpcr_all12.sh','./rmsd_sh/','./bioalign_gprot_cons_gpcr_all12/',6,50)
split_file('bioalign_gprot_all_gpcr_all12.sh','./rmsd_sh/','./bioalign_gprot_all_gpcr_all12/',6,50)
split_file('bioalign_gprot_cons_gpcr_all.sh','./rmsd_sh/','./bioalign_gprot_cons_gpcr_all/',6,50)
split_file('bioalign_gprot_all_gpcr_all.sh','./rmsd_sh/','./bioalign_gprot_all_gpcr_all/',6,50)

os.system('chmod 777 ./bioalign_gprot_cons_gpcr_all12/*.sh')
os.system('chmod 777 ./bioalign_gprot_all_gpcr_all12/*.sh')
os.system('chmod 777 ./bioalign_gprot_cons_gpcr_all/*.sh')
os.system('chmod 777 ./bioalign_gprot_all_gpcr_all/*.sh')









