#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/7_DEG/15pheatmap_GO_Astro.R --no-save > HCA_ON/scripts/7_DEG/15pheatmap_GO_Astro.out 2> HCA_ON/scripts/7_DEG/15pheatmap_GO_Astro.err
