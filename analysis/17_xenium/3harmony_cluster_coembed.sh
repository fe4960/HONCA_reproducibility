#!/bin/sh
source ~/.condainit
conda activate spatial

R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed.err 
