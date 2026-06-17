#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/13plot_impute_gene.R --no-save > HCA_ON/scripts/17_xenium/13plot_impute_gene.out 2> HCA_ON/scripts/17_xenium/13plot_impute_gene.err
