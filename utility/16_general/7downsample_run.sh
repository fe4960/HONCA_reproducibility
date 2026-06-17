#!/bin/sh

comd=HCA_ON/scripts/16_general/7downsample

mkdir $comd

fn="/dfs3b/ruic20_lab/jinjingj/human_atlas/sn/integration/sn_final_h5ad/scvi/sn_TMCB_integration_rawcount_v2"
cluster_key="subtype"
target_cells=2000

out="/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/6_downsample/sn_TMCB_integration_rawcount_v2_downsample2k_subtype"

label="majorclass"
celltype="major"

sh sb.sh -c ${comd}.sh -m 150 -t 0-1 -p fg -- $fn $cluster_key $target_cells $out $label $celltype
