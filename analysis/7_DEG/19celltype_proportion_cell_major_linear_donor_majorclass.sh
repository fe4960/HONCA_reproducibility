#!/bin/sh

comd=HCA_ON/scripts/7_DEG/19celltype_proportion_cell_major_linear_donor_majorclass

module load R/4.4.2

R < ${comd}.R --no-save > ${comd}.err 2>  ${comd}.out
