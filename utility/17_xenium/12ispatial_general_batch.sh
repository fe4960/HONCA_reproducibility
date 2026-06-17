#!/bin/sh
comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general_rest2

#comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general
mkdir $comd
for s in PP_1 #ONH_2 ONH_3 # PP_1 PP_2
do
echo $s
#sh sb.sh -c ${comd}.sh -e $comd -j $s -p fg -m 150 -t 0-20 -- $s
done


comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general_rest1
for s in PP_2
do
echo $s
#sh sb.sh -c ${comd}.sh -e $comd -j $s -p fg -m 150 -t 0-20 -- $s
done


comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/17_xenium/12ispatial_general
#for s in ONH_3 ONH_2 ONH_1 PP_1 PP_2
bc="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_xenium_w1"
for s in {1..5}
do
sh sb.sh -c ${comd}.sh -e $comd -j $s -p fg -m 150 -t 0-20 -- $s $bc
done
