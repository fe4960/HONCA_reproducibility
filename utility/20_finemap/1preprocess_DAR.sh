#!/bin/sh

#module load R/4.4.2

source ~/.condainit

conda activate archr

comd=HCA_ON/scripts/20_finemap/1preprocess_DAR

R < ${comd}.R --no-save > ${comd}.err 2> ${comd}.out
