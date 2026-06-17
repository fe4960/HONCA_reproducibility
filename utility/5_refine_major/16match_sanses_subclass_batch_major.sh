#!/bin/sh
comd="HCA_ON/scripts/5_refine_major/16match_sanses_subclass.sh"
err="HCA_ON/scripts/5_refine_major/16match_sanses_subclass"
mkdir $err

#for ct in Fibroblast BC Oligodendrocyte_precursor_cell Mural_cell AC RPE HC Cone NK_T_cell
#do

ct=major
bname=${ct}_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk_clean
sbatch -p free --mem=70GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct
#done


