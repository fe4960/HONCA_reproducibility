#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/3extract_cell_anno_bc_update_new.R --no-save > HCA_ON/scripts/17_xenium/3extract_cell_anno_bc_update_new.out 2> HCA_ON/scripts/17_xenium/3extract_cell_anno_bc_update_new.err 
