#!/bin/sh
#module load R/4.4.2 

source ~/.condainit
#conda activate spatial1
conda activate archr

R < HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap_enrichr_major.R --no-save > HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap_enrichr_major.out 2> HCA_ON/scripts/7_DEG/5GO_dotplot_comb_complexheatmap_enrichr_major.err

#R < HCA_ON/scripts/7_DEG/5GO_dotplot_comb.R --no-save > HCA_ON/scripts/7_DEG/5GO_dotplot_comb.out 2> HCA_ON/scripts/7_DEG/5GO_dotplot_comb.err
