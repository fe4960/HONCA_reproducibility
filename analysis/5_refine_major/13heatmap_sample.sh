#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/5_refine_major/13heatmap_sample.R --no-save > HCA_ON/scripts/5_refine_major/13heatmap_sample.out 2> HCA_ON/scripts/5_refine_major/13heatmap_sample.err
