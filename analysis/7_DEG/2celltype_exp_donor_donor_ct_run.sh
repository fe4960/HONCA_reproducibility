#!/bin/sh

for cell in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage 
do
zcat HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/${cell}/BCM_22_0784.txt.gz | head -n 1 | sed -e "s/\t/\n/g" > HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/ct_list_${cell}

#zcat HCA_ON/data/7_DEG/genexp_donor_subclass_final_update/${cell}/BCM_22_0047.txt.gz | head -n 1 | sed -e "s/\t/\n/g" > HCA_ON/data/7_DEG/genexp_donor_subclass_final_update/ct_list_${cell}
done
