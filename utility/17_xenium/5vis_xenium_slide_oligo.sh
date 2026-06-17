#!/bin/sh

source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/5vis_xenium_slide_oligo.R --no-save > HCA_ON/scripts/17_xenium/5vis_xenium_slide_oligo.out 2> HCA_ON/scripts/17_xenium/5vis_xenium_slide_oligo.err
