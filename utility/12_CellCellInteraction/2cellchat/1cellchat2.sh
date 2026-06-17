#!/bin/sh
source ~/.condainit
conda activate archr
R < HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat2.R --no-save > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat2.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat2.err 
