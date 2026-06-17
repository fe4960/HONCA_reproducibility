#!/bin/sh
source ~/.condainit
conda activate archr

comd="HCA_ON/scripts/14_GWAS_gene/archr_ATAC_CRE_astro_sub"

R < ${comd}.R --no-save > ${comd}.out 2> ${comd}.err
