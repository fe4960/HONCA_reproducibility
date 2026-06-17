#!/bin/sh
comd=HCA_ON/scripts/7_DEG/22sc_wilcox_DEG_30_70_mechano
mkdir $comd


c="Astrocyte"
n="Astrocyte_subclass_new5_clean"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="Fibroblast"
n="Fibroblast_subclass_clean_new"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="Oligodendrocyte"
n="Oligodendrocyte_subclass_seurat"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n



c="Endothelial_cell"
n="Endothelial_cell_subclass_sb"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="Microglia"
n="Microglia_subclass_sb"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="Macrophage"
n="Macrophage_subclass"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="MG"
n="MG_subclass_sb"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n

c="Oligodendrocyte_precursor_cell"
n="Oligodendrocyte_precursor_cell_subclass_seurat"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


c="RGC"
n="RGC_subclass_sb"

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n


for c in Cone BC HC AC 
do
n=${c}_subclass

sh sb.sh -c ${comd}.sh -p s -m 50 -t 1-10 -e $comd -j $c -- $c $n

done
