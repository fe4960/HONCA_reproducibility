#!/bin/sh
source ~/.condainit
conda activate spatial1

R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_new.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_new.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_general_new.err 
