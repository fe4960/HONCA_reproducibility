#!/bin/sh
source ~/.condainit
conda activate escape
n=$1
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general_msigdbr
mkdir $comd
R --vanilla --args $n < ${comd}.R  > ${comd}/${2}.out 2> ${comd}/${2}.err


#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general.R $n > /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general/${2}.out 2> /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general/${2}.err
