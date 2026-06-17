#!/bin/sh

comd="HCA_ON/scripts/17_xenium/10squipy_neib_region"

mkdir $comd

i=0
for g in before_LC after_LC ON mid
do
i=$(( $i + 1))
#echo "$g $i"
#sh sb.sh -c ${comd}.sh -p fg -m 50 -t 0-2 -e ${comd} -j $i -- $g $i
done


i=5
g="all"
comd="HCA_ON/scripts/17_xenium/10squipy_neib_all"
sh sb.sh -c ${comd}.sh -p fg -m 50 -t 0-2 -e ${comd} -j $i -- $g $i

