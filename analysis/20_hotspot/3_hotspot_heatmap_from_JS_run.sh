#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
comd=${main}/HCA_ON/scripts/20_hotspot/3_hotspot_heatmap_from_JS
mkdir $comd
label="Oligodendrocyte"
dir=/dfs3b/ruic20_lab/junw42/HCA_ON/data/20_hotspot/ 

#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/"
#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label

t="age"

#label="Oligodendrocyte_permut_200"
hvg="hvg_10000"
#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/"
for label in OLIGO1 OLIGO1_RBFOX1+ OLIGO1_SVEP1+ OLIGO2_LRRC7+
do
label=${label}_permute_metacell_200
echo "$label"
sh sb.sh -c ${comd}.sh -e $comd -j $label -m 50 -p fg -t 0-10 -- $dir $label $t $hvg
done

#label="Oligodendrocyte_precursor_cell_fullcell_5kgene"
label="Oligodendrocyte_precursor_cell_permut"
hvg="hvg_5000"
#hvg="hvg_10000"
#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label $t $hvg
