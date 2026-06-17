#!/bin/sh
comd="HCA_ON/scripts/5_refine_major/16match_sanses_subclass_ggplot2.sh"
err="HCA_ON/scripts/5_refine_major/16match_sanses_subclass_ggplot2"
mkdir $err

#for ct in Fibroblast BC Oligodendrocyte_precursor_cell Mural_cell AC RPE HC Cone NK_T_cell
#do

ct=Endothelial_cell
bname=${ct}_subclass_sb_seurat_rmRPE
#########sbatch -p free --mem=50GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct
#done

ct=Astrocyte
bname=${ct}_subclass_new5_clean
#######sbatch -p free --mem=50GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct

ct=Oligodendrocyte_precursor_cell
bname=${ct}_subclass_seurat_cycling
#sbatch -p free --mem=50GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct


ct=Fibroblast
bname=${ct}_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location1
#sbatch -p free --mem=60GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct


ct=Mural_cell
bname=${ct}_subclass_rmRPE
##########sbatch -p free --mem=60GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct


ct=Oligodendrocyte
bname=${ct}_subclass_seurat_rmRPE
#############sbatch -p free --mem=60GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct

ct=Macrophage
bname=${ct}_subclass_sb_clean_rmRPE
########3sbatch -p free --mem=60GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct

ct=Microglia
bname=${ct}_subclass_sb_clean
sbatch -p free --mem=60GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname $ct


ct=major
bname=${ct}_hvg10000_epochnone_seurat_rs_1_subclass_sb_seurat_rmRPE_rmUnk_clean

##########sbatch -p free --mem=20GB --time=0-1 --output=${err}/${ct}.out --error=${err}/${ct}.out $comd $bname 


