#!/bin/sh
source ~/.condainit
conda activate archr
R < HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat.R --no-save > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1cellchat.err 
