#!/bin/sh


module load R/4.4.2

R < HCA_ON/scripts/7_DEG/16example_gene_boxplot_astro.R --no-save > HCA_ON/scripts/7_DEG/16example_gene_boxplot_astro.out 2> HCA_ON/scripts/7_DEG/16example_gene_boxplot_astro.err
