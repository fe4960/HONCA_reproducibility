#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
comd=${main}/HCA_ON/scripts/20_hotspot/2_hotspot_day


label="Oligodendrocyte"

dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat.h5ad"

label=${label}_fullcell_5kgene
sh sb.sh -c ${comd}.sh -e $comd -j $label -m 200 -p fg -t 0-10 -- $dir $label


label="Oligodendrocyte_precursor_cell"

dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat.h5ad"

label=${label}_fullcell_5kgene

sh sb.sh -c ${comd}.sh -e $comd -j $label -m 200 -p fg -t 0-10 -- $dir $label

