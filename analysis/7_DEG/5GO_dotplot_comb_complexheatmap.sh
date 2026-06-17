#!/bin/sh
#module load R/4.4.2 

source ~/.condainit
#conda activate spatial1
conda activate archr

R < HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap.R --no-save > HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap.out 2> HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap.err

#R < HCA_ON/scripts/7_DEG/5GO_dotplot_comb.R --no-save > HCA_ON/scripts/7_DEG/5GO_dotplot_comb.out 2> HCA_ON/scripts/7_DEG/5GO_dotplot_comb.err
