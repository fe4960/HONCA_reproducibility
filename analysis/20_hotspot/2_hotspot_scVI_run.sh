#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
comd=${main}/HCA_ON/scripts/20_hotspot/2_hotspot_scVI


label="Oligodendrocyte"

dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat.h5ad"

label=${label}_fullcell_5kgene
#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 200 -p fg -t 1-22 -- $dir $label 5000


label="Oligodendrocyte_precursor_cell"

dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat.h5ad"



label=${label}_fullcell_10kgene

#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 200 -p fg -t 0-10 -- $dir $label 10000


label="Astrocyte"

dir=${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_new5_clean.h5ad

sh sb.sh -c ${comd}.sh -e $comd -j $label -m 70 -p s -t 1-20 -- $dir $label 5000

