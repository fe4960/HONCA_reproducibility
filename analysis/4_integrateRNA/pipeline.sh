#!/bin/sh
source ~/.condainit
conda activate scvi
main="/dfs3b/ruic20_lab/junw42"

sh ${main}/HCA_ON/scripts/4_integrateRNA/2concate2h5ad.sh 
sh ${main}/HCA_ON/scripts/4_integrateRNA/4scVI_cluster.sh
