#!/bin/sh

source ~/.condainit
conda activate scvi
python HCA_ON/scripts/7_DEG/22sc_wilcox_DEG_30_70_mechano_example.py $1 $2 $3 $4
