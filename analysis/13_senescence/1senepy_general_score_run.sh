#!/bin/sh
err=HCA_ON/scripts/13_senescence/1senepy_general_score
mkdir $err

ct="Oligodendrocyte"
fn="Oligodendrocyte_subclass_seurat"
fc=2
sh sb.sh -c ${err}.sh -p s -m 250 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc

ct=Oligodendrocyte_precursor_cell
fn=Oligodendrocyte_precursor_cell_subclass_seurat_cycling

#sh sb.sh -c ${err}.sh -p fg -m 100 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc

ct=Astrocyte
fn=Astrocyte_subclass_new5_clean

sh sb.sh -c ${err}.sh -p s -m 250 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc

ct=Microglia
fn=Microglia_subclass_sb_clean

#sh sb.sh -c ${err}.sh -p fg -m 100 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc

ct=Macrophage
fn=Macrophage_subclass_sb_clean_rmRPE

#sh sb.sh -c ${err}.sh -p fg -m 100 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc


ct=Endothelial_cell
fn=Endothelial_cell_subclass_rmRPE

#sh sb.sh -c ${err}.sh -p fg -m 100 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc

ct=Fibroblast
fn=Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location

#sh sb.sh -c ${err}.sh -p fg -m 150 -t 0-10 -e ${err} -j $ct -- $ct $fn $fc
