#!/bin/sh

comd="HCA_ON/scripts/17_xenium/7RCTD_mk"

mkdir $comd

label=oligo

cut=20
fc=5
sc=subclass
#idx="015"
idx="01"
#sh sb.sh -c ${comd}.sh -p fg -m 50 -t 0-2 -e $comd -j ${label}_${cut}_${fc}_${sc}_${idx} -- $label $cut $fc $sc $idx



cut=20
fc=5
#sc=subclass1
sc=subclass3
#idx="1" #025"
idx="n"
for label in astro #fibro astro #astro fibro
do
sh sb.sh -c ${comd}.sh -p s -m 50 -t 0-2 -e $comd -j ${label}_${cut}_${fc}_${sc}_${idx} -- $label $cut $fc $sc $idx
done
