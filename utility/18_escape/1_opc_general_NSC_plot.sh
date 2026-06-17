#!/bin/sh

source ~/.condainit
conda activate escape
#module load module load R/4.4.2
#export SSL_CERT_FILE=/etc/ssl/certs/ca-bundle.crt
#export CURL_CA_BUNDLE=/etc/ssl/certs/ca-bundle.crt
#export PATH=/pub/junw42/local/miniconda/envs/escape/bin:$PATH
#export Rscript=/pub/junw42/local/miniconda/envs/escape/bin/Rscript
n=$1
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general_NSC_plot
mkdir $comd
R --vanilla --args $n < ${comd}.R  > ${comd}/${2}.out 2> ${comd}/${2}.err


#Rscript --vanilla /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general.R $n > /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general/${2}.out 2> /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/18_escape/1_opc_general/${2}.err
