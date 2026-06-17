#!/bin/sh


comd=HCA_ON/scripts/7_DEG/1gene_count_per_donor_atlas_general
mkdir $comd

od=genexp_donor_subclass_final

ct=Astrocyte
fn=Astrocyte_subclass_new4_clean

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn  -- $ct $fn $od


ct=Oligodendrocyte
fn=${ct}_subclass_seurat

sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od

ct=Oligodendrocyte_precursor_cell
fn=${ct}_subclass_seurat

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od 
