#!/bin/sh
comd="HCA_ON/scripts/20_hotspot/2_donor_100p_metacell_permut"
mkdir $comd

label="donor"
main="/dfs3b/ruic20_lab/junw42"
cell="Oligodendrocyte"
n=200
h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_seurat"


#sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n

cell="Microglia"

#cell="Oligodendrocyte_precursor_cell"

#h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_seurat"
h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_sb_clean"


sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n

cell="Fibroblast"

h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"

sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n


cell="Endothelial_cell"


h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_rmRPE"

sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n


cell="Mural_cell"


h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_rmRPE"

sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n
