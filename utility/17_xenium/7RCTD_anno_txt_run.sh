#!/bin/sh

comd=HCA_ON/scripts/17_xenium/7RCTD_anno_txt

label="oligo"
clu="subclass"
#idx="015"
idx="01"
sh sb.sh -c ${comd}.sh -t 0-2 -p fg -m 20 -e $comd -j ${label}_${clu}_${idx} -- $label $clu $idx

label="astro"
clu="subclass1"
idx="n"

sh sb.sh -c ${comd}.sh -t 0-2 -p fg -m 20 -e $comd -j ${label}_${clu}_${idx} -- $label $clu $idx


label="fibro"
clu="subclass1"
idx="n"

sh sb.sh -c ${comd}.sh -t 0-2 -p fg -m 20 -e $comd -j ${label}_${clu}_${idx} -- $label $clu $idx

