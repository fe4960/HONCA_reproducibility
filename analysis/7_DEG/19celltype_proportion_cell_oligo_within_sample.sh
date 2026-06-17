#!/bin/sh

module load R/4.4.2

R < HCA_ON/scripts/7_DEG/19celltype_proportion_cell_oligo_within_sample.R --no-save > HCA_ON/scripts/7_DEG/19celltype_proportion_cell_oligo_within_sample.out 2> HCA_ON/scripts/7_DEG/19celltype_proportion_cell_oligo_within_sample.err
