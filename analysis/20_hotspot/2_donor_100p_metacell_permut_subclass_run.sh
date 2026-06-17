#!/bin/sh
comd="HCA_ON/scripts/20_hotspot/2_donor_100p_metacell_permut_subclass"
mkdir $comd

label="donor"
main="/dfs3b/ruic20_lab/junw42"
cell="Oligodendrocyte"
n=200
m="subclass"
h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_seurat"


sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label $n $m


cell="Oligodendrocyte_precursor_cell"

h5ad="${main}/HCA_ON/data/5_refine_major/scvi/${cell}/clean/${cell}_subclass_seurat"


#sh sb.sh -c ${comd}.sh -p fg -m 100 -t 0-5 -e $comd -j $cell -- $h5ad $cell $label
