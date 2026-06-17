#!/bin/sh


comd=HCA_ON/scripts/7_DEG/1gene_count_per_donor_atlas_general
mkdir $comd

donor="sampleid"

od=genexp_${donor}_subclass_final_new

##od=genexp_donor_subclass_final_new
sb="subclass"
ct=Astrocyte
fn=Astrocyte_subclass_new5_clean

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn  -- $ct $fn $od $sb $donor


ct=Oligodendrocyte
fn=${ct}_subclass_seurat

sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $sb $donor

ct=Oligodendrocyte_precursor_cell
fn=${ct}_subclass_seurat_cycling

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $sb $donor


ct=Fibroblast
fn=${ct}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location

#sh sb.sh -c ${comd}.sh -p fg -m 70 -t 0-2 -e $comd -j $fn -- $ct $fn $od $sb $donor
