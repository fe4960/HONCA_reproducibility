#!/bin/sh

#comd="HCA_ON/scripts/17_xenium/10squipy_neib_region"

comd="HCA_ON/scripts/17_xenium/10squipy_neib_region_ON_mod_major"

mkdir $comd

i=0
for g in before_LC after_LC ON mid
do
i=$(( $i + 1))
#echo "$g $i"
sh sb.sh -c ${comd}.sh -p fg -m 50 -t 0-2 -e ${comd} -j $i -- $g $i
done


i=3
g="ON"

#sh sb.sh -c ${comd}.sh -p fg -m 50 -t 0-2 -e ${comd} -j $i -- $g $i

