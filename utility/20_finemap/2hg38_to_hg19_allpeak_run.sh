#!/bin/sh

comd=HCA_ON/scripts/20_finemap/2hg38_to_hg19_allpeak

mkdir $comd


dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/DAR/"

list="HCA_ON/data/12_snATAC/majorclass/DAR/cell_list"

sh sb.sh -c ${comd}.sh -e $comd -j major -p s -m 10 -t 0-1 -- $dir $list
