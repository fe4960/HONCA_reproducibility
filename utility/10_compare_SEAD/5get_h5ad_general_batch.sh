#!/bin/sh

comd=HCA_ON/scripts/10_compare_SEAD/5get_h5ad_general
mkdir $comd
dir0="/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/cell_class_leiden/GSE267301"
dir1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/9_Brain_h5ad/cell_class_leiden/GSE167494"

cell=Astrocyte
fn0=${dir0}_${cell}
fn1=${dir1}_${cell}
subclass="${cell}_subclass_new4_clean"

sbatch -p free-gpu --mem=200GB --time=0-1 --account=ruic20_lab --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh $cell $fn0 $fn1 $subclass

cell=Oligodendrocyte
fn0=${dir0}_${cell}
fn1=${dir1}_${cell}
subclass="${cell}_subclass_seurat"

sbatch -p free-gpu   --mem=200GB --time=0-1 --account=ruic20_lab --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh $cell $fn0 $fn1 $subclass


cell=Oligodendrocyte_precursor_cell
fn0=${dir0}_${cell}
fn1=${dir1}_${cell}
subclass="${cell}_subclass_seurat"

sbatch -p free-gpu --mem=200GB --time=0-1 --account=ruic20_lab --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh $cell $fn0 $fn1 $subclass


cell=Microglia
fn0=${dir0}_${cell}
fn1=${dir1}_${cell}
subclass="${cell}_subclass_sb"

#sbatch --mem=200GB --time=0-1 --account=ruic20_lab --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh $cell $fn0 $fn1 $subclass

