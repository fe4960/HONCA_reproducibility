#!/bin/sh
module load R/4.4.2
R < HCA_ON/scripts/5_refine_major/12barplot_new.R --no-save > HCA_ON/scripts/5_refine_major/12barplot_new.out 2> HCA_ON/scripts/5_refine_major/12barplot_new.err
