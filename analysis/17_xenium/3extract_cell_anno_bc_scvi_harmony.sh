#!/bin/sh

source ~/.condainit
conda activate scvi

python HCA_ON/scripts/17_xenium/3extract_cell_anno_bc_scvi_harmony.py
