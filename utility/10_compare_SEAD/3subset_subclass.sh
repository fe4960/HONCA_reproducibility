#!/bin/sh
source ~/.condainit
conda activate scvi
python HCA_ON/scripts/10_compare_SEAD/3subset_subclass.py $1 $2 
