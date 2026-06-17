#!/bin/sh
source ~/.condainit
conda activate scvi
h5ad=$1
python HCA_ON/scripts/16_general/1save_obs.py $h5ad
