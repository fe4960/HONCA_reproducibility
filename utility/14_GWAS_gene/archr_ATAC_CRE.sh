#!/bin/sh
source ~/.condainit
conda activate archr

R < HCA_ON/scripts/14_GWAS_gene/archr_ATAC_CRE.R --no-save > HCA_ON/scripts/14_GWAS_gene/archr_ATAC_CRE.out 2> HCA_ON/scripts/14_GWAS_gene/archr_ATAC_CRE.err
