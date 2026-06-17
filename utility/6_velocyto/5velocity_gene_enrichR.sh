#!/bin/sh
source ~/.condainit
conda activate qtl
R < HCA_ON/scripts/6_velocyto/5velocity_gene_enrichR.R --no-save > HCA_ON/scripts/6_velocyto/5velocity_gene_enrichR.out 2> HCA_ON/scripts/6_velocyto/5velocity_gene_enrichR.err
