#!/bin/sh
source ~/.condainit
conda activate archr
main="/dfs3b/ruic20_lab/junw42"
s=$1
d=$2
Rscript --vanilla ${main}/HCA_ON/scripts/16_general/3seurat2h5ad_list.R $s $d 

