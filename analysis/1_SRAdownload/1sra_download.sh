#!/bin/sh
id=$1
output=$2

fastqdump="/pub/junw42/software/sratoolkit.3.1.1-centos_linux64/bin/fastq-dump"

#${fastqdump} --readids --split-files --gzip  -O $output $id
${fastqdump}  --split-files --gzip  -O $output $id

