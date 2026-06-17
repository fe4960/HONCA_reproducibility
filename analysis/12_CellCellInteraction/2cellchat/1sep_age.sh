#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/12_CellCellInteraction/2cellchat/1sep_age.R --no-save > HCA_ON/scripts/12_CellCellInteraction/2cellchat/1sep_age.out 2> HCA_ON/scripts/12_CellCellInteraction/2cellchat/1sep_age.err
