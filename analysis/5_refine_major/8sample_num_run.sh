#!/bin/sh
comd=HCA_ON/scripts/5_refine_major/8sample_num
mkdir $comd
s="Oligodendrocyte_precursor_cell"
fn="${s}_subclass_seurat"
#sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


s="Oligodendrocyte"
fn="${s}_subclass_seurat"
#sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn

s="Astrocyte"
#fn="${s}_subclass_new4_clean"
fn="${s}_subclass_new5_clean"

#sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


s="Fibroblast"

fn=${s}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location

sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


s="Mural_cell"

fn=${s}_subclass_rmRPE

sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


s="Macrophage"

fn=${s}_subclass_sb_clean_rmRPE

sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


s="Microglia"

fn=${s}_subclass_sb_clean

sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn


#s="Fibroblast"
#fn="${s}_subclass_clean_new"
#sh sb.sh -c ${comd}.sh -e $comd -j $s -m 10 -p fg -t 0-1 -- $s $fn
