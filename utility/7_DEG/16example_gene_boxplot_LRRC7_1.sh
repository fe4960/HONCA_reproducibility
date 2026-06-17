#!/bin/sh


module load R/4.4.2

R < HCA_ON/scripts/7_DEG/16example_gene_boxplot_LRRC7_1.R --no-save > HCA_ON/scripts/7_DEG/16example_gene_boxplot_LRRC7_1.out 2> HCA_ON/scripts/7_DEG/16example_gene_boxplot_LRRC7_1.err
