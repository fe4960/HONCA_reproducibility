#!/bin/sh
source ~/.condainit
conda activate archr
#wd=$1
#wd1=$2
#g=$3
Rscript < HCA_ON/scripts/12_CellCellInteraction/2cellchat/3build_cellchat_database.R  --no-save > HCA_ON/scripts/12_CellCellInteraction/2cellchat/3build_cellchat_database.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/3build_cellchat_database.err

#Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.R $wd $cc $ct > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.err 
