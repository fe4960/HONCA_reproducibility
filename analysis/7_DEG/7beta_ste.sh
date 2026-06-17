#!/bin/sh

for c in all #Oligodendrocyte   #Fibroblast   Astrocyte  Endothelial_cell  #all Astrocyte  Endothelial_cell  Fibroblast
do
file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${c}

perl HCA_ON/scripts/7_DEG/7beta_ste.pl $c $file
done

