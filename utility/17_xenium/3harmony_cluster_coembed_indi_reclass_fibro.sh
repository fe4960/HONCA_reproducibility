#!/bin/sh
module load R/4.4.2
R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_fibro.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_fibro.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_fibro.err
