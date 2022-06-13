#!/bin/bash
ls ./rmsd_sh/bioalign_gprot_all_gpcr_all_{1..12}.sh|xargs -n 1 -P 0 bash &&
mv *fit.txt ./bioalign_gprot_all_gpcr_all/ 
mv *.pdb ./bioalign_gprot_all_gpcr_all/ 
#ls ./rmsd_sh/bioalign_gprot_cons_gpcr_all_{1..12}.sh|xargs -n 1 -P 12 bash &&
#mv *fit.txt ./bioalign_gprot_cons_gpcr_all/ 
#mv *.pdb ./bioalign_gprot_cons_gpcr_all/
#ls ./rmsd_sh/bioalign_gprot_cons_gpcr_all12_{1..12}.sh|xargs -n 1 -P 12 bash &&
#mv *fit.txt ./bioalign_gprot_cons_gpcr_all12/
#mv *.pdb ./bioalign_gprot_cons_gpcr_all12/  
ls ./rmsd_sh/bioalign_gprot_all_gpcr_all12_{1..12}.sh|xargs -n 1 -P 12 bash &&
mv *fit.txt ./bioalign_gprot_all_gpcr_all12/ 
mv *.pdb ./bioalign_gprot_all_gpcr_all12/ 
./extract_rmsd.py  ./bioalign_gprot_all_gpcr_all/ > rmsd_gprot_all_gpcr_all.txt  
#./extract_rmsd.py  ./bioalign_gprot_cons_gpcr_all/ > rmsd_gprot_cons_gpcr_all.txt 
#./extract_rmsd.py ./bioalign_gprot_cons_gpcr_all12/ > rmsd_gprot_cons_gpcr_all12.txt 
./extract_rmsd.py ./bioalign_gprot_all_gpcr_all12/ > rmsd_gprot_all_gpcr_all12.txt 
