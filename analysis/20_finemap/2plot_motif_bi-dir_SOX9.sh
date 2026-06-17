#!/bin/sh

module load R/4.4.2
R < HCA_ON/scripts/20_finemap/2plot_motif_bi-dir_SOX9.R --no-save > HCA_ON/scripts/20_finemap/2plot_motif_bi-dir_SOX9.out 2> HCA_ON/scripts/20_finemap/2plot_motif_bi-dir_SOX9.err 

