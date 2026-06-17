#!/bin/sh
module load R/4.4.2
R < HCA_ON/scripts/14_GWAS_gene/21group_peak_bed.R --no-save > HCA_ON/scripts/14_GWAS_gene/21group_peak_bed.out 2> HCA_ON/scripts/14_GWAS_gene/21group_peak_bed.err 
