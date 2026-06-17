#!/bin/sh
#for c in all Endothelial_cell Fibroblast Oligodendrocyte Astrocyte
#for c in Oligodendrocyte 
#for c in Astrocyte
od1="genexp_donor_subclass_final_new"
label="Oligodendrocyte"
od="Oligodendrocyte_subclass_seurat"
#c=${label}

err=HCA_ON/scripts/15_manuscript/6GO_heatmap_gsea_general_GO_BP

for c in Oligodendrocyte Oligodendrocyte_precursor_cell Astrocyte Fibroblast
do
label=$c
od="${label}_subclass_seurat"
gn=300
sh sb.sh -c ${err}.sh -e ${err} -j ${label}_${gn} --  $c $od1 15 $od $label $gn
done

#label="Oligodendrocyte_precursor_cell"
#od="Oligodendrocyte_precursor_cell_subclass_seurat"
#label="Astrocyte"
#od="${label}_subclass_seurat"

###label="Fibroblast"
#####od="${label}_subclass_seurat"

######3c=${label}
######sh HCA_ON/scripts/15_manuscript/6GO_heatmap_gsea_general_GO_BP.sh $c $od1 15 $od $label
