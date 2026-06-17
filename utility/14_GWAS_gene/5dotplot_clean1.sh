#!/bin/sh
module load R/4.4.2
dir="HCA_ON/scripts/14_GWAS_gene"

s=$1
Rscript --vanilla ${dir}/5dotplot_clean1.R $s > ${dir}/5dotplot_clean1.out 2> ${dir}/5dotplot_clean1.err

#--no-save > ${dir}/5dotplot_clean1.out 2> ${dir}/5dotplot_clean1.err
