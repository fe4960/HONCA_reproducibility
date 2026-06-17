#!/bin/sh

source ~/.condainit

conda activate scvi

python HCA_ON/scripts/5_refine_major/15generate_majorclass.py
