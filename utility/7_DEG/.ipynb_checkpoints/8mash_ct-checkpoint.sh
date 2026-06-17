#!/bin/sh
source ~/.condainit
conda activate qtl
R < HCA_ON/scripts/7_DEG/8mash_ct.R --no-save > HCA_ON/scripts/7_DEG/8mash_ct.out 2> HCA_ON/scripts/7_DEG/8mash_ct.err 
