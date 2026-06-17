#!/bin/sh
source ~/.condainit
conda activate archr
wd=$1
cc=$2
ct=$3
h5ad=$4
Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.R $wd $cc $ct $h5ad

#Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.R $wd $cc $ct > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.err 
