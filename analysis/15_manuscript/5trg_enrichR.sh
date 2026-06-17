#!/bin/sh
source ~/.condainit
conda activate qtl
#comd="HCA_ON/scripts/15_manuscript/5trg_enrichR_opc_cyc"

comd="HCA_ON/scripts/15_manuscript/5trg_enrichR"
R < ${comd}.R --no-save > ${comd}.out 2> ${comd}.err
