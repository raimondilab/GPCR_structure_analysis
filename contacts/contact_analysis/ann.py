
import os
import pickle
#####Annotate contact interface for each structure  with GPCRDB numbering and CGN numbering
# input data
infile = open('./use_file/Gprotein_position.pickle', 'rb')
gprotein = pickle.load(infile)
infile1 = open('./use_file/GPCRDB_pos.pickle', 'rb')
gpc = pickle.load(infile1)
input_files = [f for f in os.listdir('./cont_file') if f.endswith('.tsv')]

for input_file in input_files:
    pdb = input_file.split("_")[0]
    gpcr1 = input_file.split("_")[1]
    gprot2 = input_file.split("_")[2]
    output = './ann_file/' + pdb + "_" + gpcr1 + "_" + gprot2 + "_2.txt"
    
    # mapping
    t = {}
    with open('./cont_file/'+input_file) as f:
        next(f)
        for line in f:
            gpcr, gpcr_uniprot, gpcr_pos, gprot, gprot_id, gprot_pos = line.strip().split('\t')
            gprot_info = gprotein.get((gprot_id, gprot_pos))
            if gprot_info is None:
                gprot_info = ["", ""]
            gpc_info = gpc.get((gpcr_uniprot, int(gpcr_pos)))
            if gpc_info is None:
                gpc_info = ["", "", ""]
            t[(gpcr, gpcr_pos, gprot, gprot_pos)] = [gpcr_uniprot, gprot_info[0], gprot_info[1], gpc_info[1], gpc_info[2]]
    
    # save output
    with open(output, "w") as f:
        print("\t".join(['GPCR', 'Uniprot', 'Pos1', 'BW', 'Structural', 'Gprotein', 'Pos2', 'Gprot_pos', 'Gprot_struct']), file=f)
        for k, v in t.items():
            print("\t".join([k[0], v[0], k[1], v[3], v[4], k[2], k[3], v[1], v[2]]), file=f)
