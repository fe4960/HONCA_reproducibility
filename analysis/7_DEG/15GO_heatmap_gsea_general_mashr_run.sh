#!/bin/sh
for c in all #Endothelial_cell Fibroblast Oligodendrocyte Astrocyte
do
sh HCA_ON/scripts/7_DEG/15GO_heatmap_gsea_general_mashr.sh $c age 10
done
