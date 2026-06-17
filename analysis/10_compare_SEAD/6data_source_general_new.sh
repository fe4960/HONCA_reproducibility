#!/bin/sh
source ~/.condainit
conda activate scvi
cell=$1
python HCA_ON/scripts/10_compare_SEAD/6data_source_general_new.py $cell
