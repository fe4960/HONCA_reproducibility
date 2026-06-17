#!/bin/sh
source ~/.condainit
conda activate archr
s=$1
fn=$2
Rscript --vanilla HCA_ON/scripts/5_refine_major/8sample_num.R $s $fn > HCA_ON/scripts/5_refine_major/8sample_num.out 2> HCA_ON/scripts/5_refine_major/8sample_num.err
