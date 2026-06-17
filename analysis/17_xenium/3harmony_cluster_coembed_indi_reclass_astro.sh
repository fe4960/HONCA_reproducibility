#!/bin/sh
module load R/4.4.2
R < HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_astro.R --no-save > HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_astro.out 2> HCA_ON/scripts/17_xenium/3harmony_cluster_coembed_indi_reclass_astro.err
