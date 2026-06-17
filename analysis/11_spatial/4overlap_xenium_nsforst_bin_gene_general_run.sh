#!/bin/sh

#label="subclass"
label="majorclass"
pa="hecaXenium_478.txt"
#pa="Xenium_hLung_v1_metadata.txt"
for ct in Fibroblast #Oligodendrocyte Oligodendrocyte_precursor_cell Astrocyte
do
sh HCA_ON/scripts/11_spatial/4overlap_xenium_nsforst_bin_gene_general.sh $ct $label $pa
done
