#!/bin/sh
source ~/.condainit
conda activate archr
wd=$1
wd1=$2
g=$3
R --vanilla --args  $wd ${wd1} $g < HCA_ON/scripts/12_CellCellInteraction/2cellchat/4plot_astro.R 

#Rscript --vanilla HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.R $wd $cc $ct > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat1_general.err 
