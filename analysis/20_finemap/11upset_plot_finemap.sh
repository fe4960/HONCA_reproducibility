#!/bin/sh

source ~/.condainit

conda activate archr

R < HCA_ON/scripts/20_finemap/11upset_plot_finemap.R --no-save > HCA_ON/scripts/20_finemap/11upset_plot_finemap.out 2> HCA_ON/scripts/20_finemap/11upset_plot_finemap.err 
