#!/bin/sh

source ~/.condainit
conda activate scvi

python HCA_ON/scripts/5_refine_major/10obs_sum.py
