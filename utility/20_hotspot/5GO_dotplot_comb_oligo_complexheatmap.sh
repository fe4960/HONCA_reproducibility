#!/bin/sh
#module load R/4.4.2 

source ~/.condainit
conda activate archr

R < HCA_ON/scripts/20_hotspot/5GO_dotplot_comb_oligo_complexheatmap.R --no-save > HCA_ON/scripts/20_hotspot/5GO_dotplot_comb_oligo_complexheatmap.out 2> HCA_ON/scripts/20_hotspot/5GO_dotplot_comb_oligo_complexheatmap.err
