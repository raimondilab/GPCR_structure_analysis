
import pandas as pd
import os,sys
import numpy as np
import pickle
##########Calculate contact fractions for each Gprotein family and network link(edge) and node files

df=pd.read_csv("./use_file/all_gpcrs.tsv", comment="#", sep="\t")

equiv={'GNAS':'Gs','GNAL':'Gs','GNAI1':'Gio','GNAI2':'Gio','GNAI3':'Gio','GNAO1':'Gio','GNAZ':'Gio','GNAT1':'Gio','GNAT2':'Gio','GNA15':'Gq11','GNA11':'Gq11','GNAQ':'Gq11','GNA14':'Gq11','GNA12':'G1213','GNA13':'G1213'}
df["gfam"] = df["Gprotein"].map(equiv)
df['pair']=df["BW"]+"|"+df["Gprot_pos"]
df['pair_struct']=df["Structural"] +"|"+ df["Gprot_struct"]
gs=df[df.gfam =='Gs']
gio=df[df.gfam =='Gio']
gq11=df[df.gfam =='Gq11']
g1213=df[df.gfam =='G1213']

#gs
gs['unique_gpcr_contact'] = gs.groupby(['pair'])['GPCR'].transform('nunique')
gs['unique_gpcr_contact_tm'] = gs.groupby(['pair_struct'])['GPCR'].transform('nunique')
gs['fraction_pair']=gs['unique_gpcr_contact']/len(set(gs['GPCR']))
#gs = gs[(gs['unique_gpcr_contact']/len(set(gs['GPCR'])) > 0.25)]rm('mean')
gs['fraction_pair_struct']=gs.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gs['num_cont_gprot'] = gs.groupby(['Gprot_struct'])['pair'].transform('nunique')
gs['num_cont_gpcr'] = gs.groupby(['Structural'])['pair'].transform('nunique') 

gs.to_csv('./network/files/gs_contact_all.tsv',index=None,header=True,sep='\t')
#gio
gio['unique_gpcr_contact'] = gio.groupby(['pair'])['GPCR'].transform('nunique')
gio['unique_gpcr_contact_tm'] = gio.groupby(['pair_struct'])['GPCR'].transform('nunique')
gio['fraction_pair']=gio['unique_gpcr_contact']/len(set(gio['GPCR']))
#gio = gio[(gio['unique_gpcr_contact']/len(set(gio['GPCR'])) > 0.25)]
gio['fraction_pair_struct']=gio.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gio['num_cont_gprot'] = gio.groupby(['Gprot_struct'])['pair'].transform('nunique')
gio['num_cont_gpcr'] = gio.groupby(['Structural'])['pair'].transform('nunique') 
gio.to_csv('./network/files/gio_contact_all.tsv',index=None,header=True,sep='\t')

#gs
nodes1=gs[["Structural","num_cont_gpcr"]].drop_duplicates()
nodes2=gs[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
nodes1.columns = ['struct','num']
nodes2.columns = ['struct','num']
nodes = pd.concat([nodes1, nodes2], axis=0)
nodes.columns = ['struct','num_cont_gs']
nodes.to_csv('./network/files/gs_nodes_all.tsv',index=None,header=True,sep='\t')
links=gs[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
links.columns = ['SS1','SS2','cont_gs',"fraction_gs"]
links.to_csv('./network/files/gs_links_all.tsv',index=None,header=True,sep='\t')

#gio
nodes4=gio[["Structural","num_cont_gpcr"]].drop_duplicates()
nodes5=gio[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
nodes4.columns = ['struct','num']
nodes5.columns = ['struct','num']
nodes3 = pd.concat([nodes4, nodes5], axis=0)
nodes3.columns = ['struct','num_cont_gio']
nodes3.to_csv('./network/files/gio_nodes_all.tsv',index=None,header=True,sep='\t')
links2=gio[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
links2.columns = ['SS1','SS2','cont_gio',"fraction_gio"]
links2.to_csv('./network/files/gio_links_all.tsv',index=None,header=True,sep='\t')

#gq11
gq11['unique_gpcr_contact'] = gq11.groupby(['pair'])['GPCR'].transform('nunique')
gq11['unique_gpcr_contact_tm'] = gq11.groupby(['pair_struct'])['GPCR'].transform('nunique')
gq11['fraction_pair']=gq11['unique_gpcr_contact']/len(set(gq11['GPCR']))
#gq11 = gq11[(gq11['unique_gpcr_contact']/len(set(gq11['GPCR'])) > 0.25)]
gq11['fraction_pair_struct']=gq11.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gq11['num_cont_gprot'] = gq11.groupby(['Gprot_struct'])['pair'].transform('nunique')
gq11['num_cont_gpcr'] = gq11.groupby(['Structural'])['pair'].transform('nunique') 
gq11.to_csv('./network/files/gq11_contact_all.tsv',index=None,header=True,sep='\t')


#g1213
g1213['unique_gpcr_contact'] = g1213.groupby(['pair'])['GPCR'].transform('nunique')
g1213['unique_gpcr_contact_tm'] = g1213.groupby(['pair_struct'])['GPCR'].transform('nunique')
g1213['fraction_pair']=g1213['unique_gpcr_contact']/len(set(g1213['GPCR']))
#g1213 = g1213[(g1213['unique_gpcr_contact']/len(set(g1213['GPCR'])) > 0.25)]
g1213['fraction_pair_struct']=g1213.groupby(['pair_struct'])['fraction_pair'].transform('mean')
g1213['num_cont_gprot'] = g1213.groupby(['Gprot_struct'])['pair'].transform('nunique')
g1213['num_cont_gpcr'] = g1213.groupby(['Structural'])['pair'].transform('nunique') 
g1213.to_csv('./network/files/g1213_contact_all.tsv',index=None,header=True,sep='\t')

#s=pd.concat([gs,gio,gq11,g1213])
s=pd.concat([gs,gio,gq11,g1213])
#s=pd.concat([gs,gio])
s=s[['BW','Gprot_pos', 'gfam','fraction_pair']].drop_duplicates()
s.to_csv('./network/files/fraction_pair_all.tsv',index=None,header=True,sep='\t')

#gq11
nod1=gq11[["Structural","num_cont_gpcr"]].drop_duplicates()
nod2=gq11[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
nod1.columns = ['struct','num']
nod2.columns = ['struct','num']
nod = pd.concat([nod1, nod2], axis=0)
nod.columns = ['struct','num_cont_gq11']
nod.to_csv('./network/files/gq11_nodes_all.tsv',index=None,header=True,sep='\t')
lin=gq11[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
lin.columns = ['SS1','SS2','cont_gq11',"fraction_gq11"]
lin.to_csv('./network/files/gq11_links_all.tsv',index=None,header=True,sep='\t')

#g1213
nod4=g1213[["Structural","num_cont_gpcr"]].drop_duplicates()
nod5=g1213[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
nod4.columns = ['struct','num']
nod5.columns = ['struct','num']
nod3 = pd.concat([nod4, nod5], axis=0)
nod3.columns = ['struct','num_cont_g1213']
nod3.to_csv('./network/files/g1213_nodes_all.tsv',index=None,header=True,sep='\t')
lin2=g1213[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
lin2.columns = ['SS1','SS2','cont_g1213',"fraction_g1213"]
lin2.to_csv('./network/files/g1213_links_all.tsv',index=None,header=True,sep='\t')


####classA
infile = open('./use_file/GPCR_class.pickle', 'rb')
clas = pickle.load(infile)
cla=pd.DataFrame(clas.items(), columns=['GPCR', 'class'])
classA=df.merge(cla,left_on=['GPCR'],right_on=['GPCR'],how="left")
classA=classA.loc[classA['class'] == 'classA']
gs1=classA[classA.gfam =='Gs']
gio1=classA[classA.gfam =='Gio']
gq111=classA[classA.gfam =='Gq11']
g12131=classA[classA.gfam =='G1213']


#gs1
gs1['unique_gpcr_contact'] = gs1.groupby(['pair'])['GPCR'].transform('nunique')
gs1['unique_gpcr_contact_tm'] = gs1.groupby(['pair_struct'])['GPCR'].transform('nunique')
gs1['fraction_pair']=gs1['unique_gpcr_contact']/len(set(gs1['GPCR']))
#gs1 = gs1[(gs1['unique_gpcr_contact']/len(set(gs1['GPCR'])) > 0.25)]
gs1['fraction_pair_struct']=gs1.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gs1['num_cont_gprot'] = gs1.groupby(['Gprot_struct'])['pair'].transform('nunique')
gs1['num_cont_gpcr'] = gs1.groupby(['Structural'])['pair'].transform('nunique') 
gs1.to_csv('./network/files/gs_contact_all_classA.tsv',index=None,header=True,sep='\t')
#gio1
gio1['unique_gpcr_contact'] = gio1.groupby(['pair'])['GPCR'].transform('nunique')
gio1['unique_gpcr_contact_tm'] = gio1.groupby(['pair_struct'])['GPCR'].transform('nunique')
gio1['fraction_pair']=gio1['unique_gpcr_contact']/len(set(gio1['GPCR']))
#gio1 = gio1[(gio1['unique_gpcr_contact']/len(set(gio1['GPCR'])) > 0.25)]
gio1['fraction_pair_struct']=gio1.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gio1['num_cont_gprot'] = gio1.groupby(['Gprot_struct'])['pair'].transform('nunique')
gio1['num_cont_gpcr'] = gio1.groupby(['Structural'])['pair'].transform('nunique') 
gio1.to_csv('./network/files/gio_contact_all_classA.tsv',index=None,header=True,sep='\t')

#gs1
node1=gs1[["Structural","num_cont_gpcr"]].drop_duplicates()
node2=gs1[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
node1.columns = ['struct','num']
node2.columns = ['struct','num']
node = pd.concat([node1, node2], axis=0)
node.columns = ['struct','num_cont_gs']
node.to_csv('./network/files/gs_nodes_all_classA.tsv',index=None,header=True,sep='\t')
link=gs1[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
link.columns = ['SS1','SS2','cont_gs',"fraction_gs"]
link.to_csv('./network/files/gs_links_all_classA.tsv',index=None,header=True,sep='\t')

#gio1
node4=gio1[["Structural","num_cont_gpcr"]].drop_duplicates()
node5=gio1[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
node4.columns = ['struct','num']
node5.columns = ['struct','num']
node3 = pd.concat([node4, node5], axis=0)
node3.columns = ['struct','num_cont_gio']
node3.to_csv('./network/files/gio_nodes_all_classA.tsv',index=None,header=True,sep='\t')
link2=gio1[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
link2.columns = ['SS1','SS2','cont_gio',"fraction_gio"]
link2.to_csv('./network/files/gio_links_all_classA.tsv',index=None,header=True,sep='\t')


#gq111
gq111['unique_gpcr_contact'] = gq111.groupby(['pair'])['GPCR'].transform('nunique')
gq111['unique_gpcr_contact_tm'] = gq111.groupby(['pair_struct'])['GPCR'].transform('nunique')
gq111['fraction_pair']=gq111['unique_gpcr_contact']/len(set(gq111['GPCR']))
#gq111 = gq111[(gq111['unique_gpcr_contact']/len(set(gq111['GPCR'])) > 0.25)]
gq111['fraction_pair_struct']=gq111.groupby(['pair_struct'])['fraction_pair'].transform('mean')
gq111['num_cont_gprot'] = gq111.groupby(['Gprot_struct'])['pair'].transform('nunique')
gq111['num_cont_gpcr'] = gq111.groupby(['Structural'])['pair'].transform('nunique') 
gq111.to_csv('./network/files/gq11_contact_all_classA.tsv',index=None,header=True,sep='\t')


#g12131
g12131['unique_gpcr_contact'] = g12131.groupby(['pair'])['GPCR'].transform('nunique')
g12131['unique_gpcr_contact_tm'] = g12131.groupby(['pair_struct'])['GPCR'].transform('nunique')
g12131['fraction_pair']=g12131['unique_gpcr_contact']/len(set(g12131['GPCR']))
#g12131 = g12131[(g12131['unique_gpcr_contact']/len(set(g12131['GPCR'])) > 0.25)]
g12131['fraction_pair_struct']=g12131.groupby(['pair_struct'])['fraction_pair'].transform('mean')
g12131['num_cont_gprot'] = g12131.groupby(['Gprot_struct'])['pair'].transform('nunique')
g12131['num_cont_gpcr'] = g12131.groupby(['Structural'])['pair'].transform('nunique') 
g12131.to_csv('./network/files/g1213_contact_all_classA.tsv',index=None,header=True,sep='\t')

s1=pd.concat([gs1,gq111,g12131,g12131])
#s1=pd.concat([gs1,gq111,gio1])
#s1=pd.concat([gs1,gio1])
s1=s1[['BW','Gprot_pos', 'gfam','fraction_pair']].drop_duplicates()
s1.to_csv('./network/files/fraction_pair__all_classA.tsv',index=None,header=True,sep='\t')

#gq111
no1=gq111[["Structural","num_cont_gpcr"]].drop_duplicates()
no2=gq111[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
no1.columns = ['struct','num']
no2.columns = ['struct','num']
no = pd.concat([no1, no2], axis=0)
no.columns = ['struct','num_cont_gq11']
no.to_csv('./network/files/gq11_nodes_all_classA.tsv',index=None,header=True,sep='\t')
li=gq111[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
li.columns = ['SS1','SS2','cont_gq11',"fraction_gq11"]
li.to_csv('./network/files/gq11_links_all_classA.tsv',index=None,header=True,sep='\t')


#g1213
no4=g12131[["Structural","num_cont_gpcr"]].drop_duplicates()
no5=g12131[["Gprot_struct","num_cont_gprot"]].drop_duplicates()
no4.columns = ['struct','num']
no5.columns = ['struct','num']
no3 = pd.concat([node4, node5], axis=0)
no3.columns = ['struct','num_cont_gio']
no3.to_csv('./network/files/g1213_nodes_all_classA.tsv',index=None,header=True,sep='\t')
li2=g12131[['Structural','Gprot_struct','unique_gpcr_contact_tm','fraction_pair_struct']].drop_duplicates()
li2.columns = ['SS1','SS2','cont_g1213',"fraction_g1213"]
li2.to_csv('./network/files/g1213_links_all_classA.tsv',index=None,header=True,sep='\t')




