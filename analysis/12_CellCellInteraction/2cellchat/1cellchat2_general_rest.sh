#!/bin/sh
source ~/.condainit
conda activate archr
wd=$1
#cc=$2
Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat2_general_rest.R $wd

#Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.R $wd $cc $ct > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.err 
