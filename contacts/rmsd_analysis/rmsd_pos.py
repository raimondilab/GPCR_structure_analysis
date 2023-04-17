
import pandas as pd
import numpy as np
import pickle

######### Extracting common positions for rmsd calculation for GPCR and G protein part of the structures

df=pd.read_csv("/home/mmatic/GPCR_structure_analysis/relax_rmsd/use_file/GPCR_structs_clean.tsv", comment="#", sep="\t")
infile = open('/home/mmatic/GPCR_structure_analysis/relax_rmsd/use_file/Gprotein_position.pickle', 'rb')
gprotein = pickle.load(infile)
infile1 = open('/home/mmatic/GPCR_structure_analysis/relax_rmsd/use_file/GPCRDB_pos.pickle', 'rb')
gpc = pickle.load(infile1)

##### mapping
gprot={}
for index,row in df.iterrows():
	gprot[(row['Gprotein_Uniprot_AC'],row['Gprotein_Uniprot_ID'])]=row['Gprotein_Gene_name']
gp={}
for k,v in gprot.items():
    for k1,v1 in gprotein.items():
        if k[1]==k1[0]:
            gp[(k[0],int(k1[1]))]=[v1[0],v1[1],v,k1[0]]
        
    
####### mapping
    t = {}
    with open('/home/mmatic/GPCR_structure_analysis/relax_rmsd/use_file/uniprot_pdb.tsv') as f:
        next(f)
        for line in f:
            pdbid, chain, pdb_pos, uniprot_id, uniprot_pos = line.strip().split('\t')
            gprot_info = gp.get((uniprot_id,  int(uniprot_pos)))
            if gprot_info is None:
                gprot_info = ["", "","",""]
            gpc_info = gpc.get((uniprot_id, int(uniprot_pos)))
            if gpc_info is None:
                gpc_info = ["", "", ""]
            t[(pdbid, chain, uniprot_id,uniprot_pos)] = [pdb_pos,gprot_info[0], gprot_info[1],gprot_info[2],gpc_info[0], gpc_info[1], gpc_info[2]]
    

#####preparing dataframe
t1=pd.DataFrame.from_dict(t, orient='index').reset_index()
t1[['b1', 'b2','b3','b4']] = pd.DataFrame(t1['index'].tolist(), index=t1.index)                                                                                                
t1=t1.drop(['index'], axis=1)
t1.columns = ['PDB_pos', 'Gprot_pos', 'Gprot_struct', 'Gprotein', 'GPCR', 'BW', 'Structural', 'PDBID', 'Chain', 'Uniprot_ID', 'Uniprot_pos']
t1 = t1.replace(r'^\s*$', np.nan, regex=True)
t1 = t1.dropna(subset=['GPCR', 'Gprotein'], how='all')
####correcting wrong chain mapping
chains=(t1['PDBID']+'-'+t1['Chain']+'-'+t1['Uniprot_ID']).unique().tolist()
chains1=(df['PDB_ID']+'-'+df['Receptor_Chain']+'-'+df['Receptor_Uniprot_AC']).unique().tolist()+(df['PDB_ID']+'-'+df['Gprotein_chain']+'-'+df['Gprotein_Uniprot_AC']).unique().tolist()
list(set(chains) - set(chains1))
chain={}
for i in chains:
    for j in chains1:
        if i.split('-')[0] in j.split('-')[0] and i.split('-')[2] in j.split('-')[2] and i.split('-')[1] not in j.split('-')[1]:
            chain[i.split('-')[0]+'-'+i.split('-')[1]]=j.split('-')[0]+'-'+j.split('-')[1]
for k ,v in chain.items():
    for index,row in t1.iterrows():
        if k.split('-')[0] in row['PDBID'] and k.split('-')[1] in row['Chain']:
            row['Chain'] = row['Chain'].replace(row['Chain'],v.split('-')[1])
t1 = t1[~((t1['PDBID'] == '7eb2') & (t1['Chain'] == 'C'))]
#####extracting positions
gpcrs = t1[['PDB_pos','BW', 'Structural', 'GPCR','PDBID', 'Chain', 'Uniprot_ID',]].dropna()
gpcrs=gpcrs[~(gpcrs['BW'] == '2.58')]
gpcrs['Chain'] = gpcrs['Chain'].str.replace('RP', 'E')

gproteins = t1[['PDB_pos','Gprot_pos', 'Gprot_struct', 'Gprotein','PDBID', 'Chain', 'Uniprot_ID',]].dropna()
gproteins['Chain'] = gproteins['Chain'].str.replace('AP', 'A')

###gpcr all
d5=gpcrs.groupby(['BW'])['PDBID'].count().reset_index()
d5=d5.loc[d5['PDBID'].isin([df['PDB_ID'].unique().size])]
d5.to_csv('./pos/gpcr_pos_set_all.tsv',index=None,header=True,sep='\t')
df6=gpcrs.loc[gpcrs['BW'].isin(['1.49','1.50','1.51','2.49','2.50','2.51','3.49','3.50','3.51','4.49','4.50','4.51','5.49','5.50','5.51','6.49','6.50','7.49','7.50','7.51'])]
##gprotein
d6=gproteins.groupby(['Gprot_pos'])["PDBID"].count().reset_index()
d6=d6.loc[d6['PDBID'].isin([df['PDB_ID'].unique().size])]
d6.to_csv('./pos/gprot_pos_set_all.tsv',index=None,header=True,sep='\t')
df7=gproteins.loc[gproteins['Gprot_pos'].isin(d6['Gprot_pos'].tolist())]
#gpcr all
df8=gpcrs.loc[gpcrs['BW'].isin(d5['BW'].tolist())]
df6['PDB_pos']= df6['PDB_pos'].astype(str)
df7['PDB_pos']= df7['PDB_pos'].astype(str)
df8['PDB_pos']= df8['PDB_pos'].astype(str)
ds6=df6.groupby(['PDBID','Chain'])["PDB_pos"].apply(lambda item:','.join(item)).reset_index()
ds7=df7.groupby(['PDBID','Chain'])["PDB_pos"].apply(lambda item:','.join(item)).reset_index()
ds8=df8.groupby(['PDBID','Chain'])["PDB_pos"].apply(lambda item:','.join(item)).reset_index()
ds6.to_csv('./pos/gpcr_pos_all_12.tsv',index=None,header=None,sep='\t')
ds7.to_csv('./pos/gprot_pos_all.tsv',index=None,header=None,sep='\t')
ds8.to_csv('./pos/gpcr_pos_all.tsv',index=None,header=None,sep='\t')
