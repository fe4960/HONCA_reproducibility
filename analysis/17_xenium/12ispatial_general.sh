#!/bin/sh

source ~/.condainit
conda activate spatial1

#comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general_astro
mkdir $comd
s=$1
b=$2
R --vanilla --args $s $b < ${comd}.R > ${comd}/${s}.out 2> ${comd}/${s}.err

#R < ${comd}.R $s > ${comd}/${s}.out 2> ${comd}/${s}.err 
