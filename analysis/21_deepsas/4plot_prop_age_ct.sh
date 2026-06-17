#!/bin/sh

module load R/4.4.2 

R < HCA_ON/scripts/21_deepsas/4plot_prop_age_ct.R  --no-save > HCA_ON/scripts/21_deepsas/4plot_prop_age_ct.out 2> HCA_ON/scripts/21_deepsas/4plot_prop_age_ct.err
