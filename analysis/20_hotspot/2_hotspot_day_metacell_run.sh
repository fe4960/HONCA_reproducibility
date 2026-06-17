#!/bin/sh
main="/dfs3b/ruic20_lab/junw42"
comd=${main}/HCA_ON/scripts/20_hotspot/2_hotspot_day_metacell


label="Oligodendrocyte"

#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat_leiden_${label}_permute_metacell.h5ad"
dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat_leiden_${label}_permute_metacell_200.h5ad"
hvg=10000
label=${label}_permut_200
###dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat_leiden_donor_metacell.h5ad"
#########sh sb.sh -c ${comd}.sh -e $comd -j $label -m 200 -p fg -t 0-10 -- $dir $label $hvg


label="Oligodendrocyte_precursor_cell"
dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat_leiden_${label}_permute_metacell.h5ad"
label=${label}_permut
hvg=10000
#dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_seurat_leiden_donor_metacell.h5ad"
########sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label $hvg

hvg=10000

for c in OLIGO1 OLIGO1_RBFOX1+ OLIGO1_SVEP1+ OLIGO2_LRRC7+
do
label=${c}_permute_metacell_200
dir="${main}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_leiden_${label}.h5ad"
#sh sb.sh -c ${comd}.sh -e $comd -j $label -m 50 -p fg -t 0-5 -- $dir $label $hvg
done


label="Microglia"
dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_sb_clean_leiden_${label}_permute_metacell_200.h5ad"
label=${label}_permut
hvg=10000

sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label $hvg

label="Fibroblast"
dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location_leiden_${label}_permute_metacell_200.h5ad"
label=${label}_permut
hvg=10000

sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label $hvg

for label in Mural_cell Endothelial_cell
do
dir="${main}/HCA_ON/data/5_refine_major/scvi/${label}/clean/${label}_subclass_rmRPE_leiden_${label}_permute_metacell_200.h5ad"
label=${label}_permut
hvg=10000
sh sb.sh -c ${comd}.sh -e $comd -j $label -m 150 -p fg -t 0-10 -- $dir $label $hvg
done
