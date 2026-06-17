#!/bin/sh
#for c in all Endothelial_cell Fibroblast Oligodendrocyte Astrocyte
#for c in Oligodendrocyte 
#for c in Astrocyte
od="genexp_donor_subclass_final"
label="Oligodendrocyte_2"
for c in major #Oligodendrocyte
do
sh HCA_ON/scripts/7_DEG/15GO_heatmap_gsea_general_major.sh $c age_cpm1 10 #$od $label

	#sh HCA_ON/scripts/7_DEG/15GO_heatmap_gsea_general.sh $c age 10
done
