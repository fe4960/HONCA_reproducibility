#!/bin/sh

source ~/.condainit

conda activate archr

R < HCA_ON/scripts/7_DEG/11upset_plot_age_nomashr_sub.R --no-save > HCA_ON/scripts/7_DEG/11upset_plot_age_nomashr_sub.out 2> HCA_ON/scripts/7_DEG/11upset_plot_age_nomashr_sub.err
