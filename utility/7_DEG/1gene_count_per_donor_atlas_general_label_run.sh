#!/bin/sh


comd=HCA_ON/scripts/7_DEG/1gene_count_per_donor_atlas_general_label
mkdir $comd

od=genexp_donor_subclass_final_update
#od=genexp_donor_major_final
lab="subclass"

ct=Astrocyte
fn=Astrocyte_subclass_new5_clean

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn  -- $ct $fn $od $lab


ct=Oligodendrocyte
fn=${ct}_subclass_seurat

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab

ct=Oligodendrocyte_precursor_cell
fn=${ct}_subclass_seurat

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab


ct="Fibroblast"
fn="Fibroblast_subclass_clean_new"

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab

ct="Endothelial_cell"
fn="Endothelial_cell_subclass_sb"

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab


ct="Microglia"
fn="Microglia_subclass_sb"

sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab


ct="Macrophage"
fn="Macrophage_subclass"

sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab

ct="MG"
fn="MG_subclass_sb"


sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab


ct="RGC"
fn="RGC_subclass_sb"


sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $lab

