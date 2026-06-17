#!/bin/sh

zcat ${1} | awk '{if (NR % 4 == 0) { a=length($0); if(a==26){n1=n1+1;}else if(a==27){n2=n2+1;}else if(a==28){n3=n3+1;} ; print n1, n2, n3}}' | tail -n 1 | awk '{total=$1+$2+$3; b1=$1/total; b2=$2/total; b3=$3/total; print "26bp:"$1"\t"b1"\n27bp:"$2"\t"b2"\n28bp:"$3"\t"b3"\n" }'  > ${1}.len
