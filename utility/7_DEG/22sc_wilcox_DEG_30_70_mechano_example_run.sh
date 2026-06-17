#!/bin/sh

ct="Oligodendrocyte"
name="Oligodendrocyte_subclass_seurat"
st="OLIGO2_LRRC7+"
gene="PIEZO2"


sh HCA_ON/scripts/7_DEG/22sc_wilcox_DEG_30_70_mechano_example.sh $ct $name $st $gene
