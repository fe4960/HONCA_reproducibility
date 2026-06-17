#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
comd=${main}/HCA_ON/scripts/20_hotspot/3_hotspot_heatmap
mkdir $comd
label="Oligodendrocyte"
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/ 

#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/"
#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label


#label="Oligodendrocyte_precursor_cell"
label="Oligodendrocyte_precursor_cell_fullcell_5kgene"

#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/"
sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label

