#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/7RCTD_major_anno_xenium_clu.R --no-save > HCA_ON/scripts/17_xenium/7RCTD_major_anno_xenium_clu.out 2> HCA_ON/scripts/17_xenium/7RCTD_major_anno_xenium_clu.err 
