#!/bin/sh
module load R/4.4.2

R < HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_opc_sample_comb.R --no-save > HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_opc_sample_comb.out 2> HCA_ON/scripts/7_DEG/19celltype_proportion_cell_sub_opc_sample_comb.err
