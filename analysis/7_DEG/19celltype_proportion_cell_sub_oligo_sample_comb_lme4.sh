#!/bin/sh
module load R/4.4.2

R < HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_oligo_sample_comb_lme4.R --no-save > HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_oligo_sample_comb_lme4.out 2> HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_oligo_sample_comb_lme4.err
