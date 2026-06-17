#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/3extract_cell_anno_bc.R --no-save > HCA_ON/scripts/17_xenium/3extract_cell_anno_bc.out 2> HCA_ON/scripts/17_xenium/3extract_cell_anno_bc.err 
