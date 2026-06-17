#!/bin/sh
source ~/.condainit
conda activate qtl

R < HCA_ON/scripts/14_GWAS_gene/5dotplot_clean_LDSC.R --no-save > HCA_ON/scripts/14_GWAS_gene/5dotplot_clean_LDSC.out 2> HCA_ON/scripts/14_GWAS_gene/5dotplot_clean_LDSC.err
