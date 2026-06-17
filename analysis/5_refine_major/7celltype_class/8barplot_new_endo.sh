#!/bin/sh
module load R/4.2.2

R < HCA_ON/scripts/5_refine_major/7celltype_class/8barplot_new_endo.R --no-save > HCA_ON/scripts/5_refine_major/7celltype_class/8barplot_new_endo.out 2> HCA_ON/scripts/5_refine_major/7celltype_class/8barplot_new_endo.err
