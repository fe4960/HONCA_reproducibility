#!/bin/sh

comd=HCA_ON/scripts/16_general/7downsample

mkdir $comd

#fn="/dfs3b/ruic20_lab/jinjingj/human_atlas/sn/integration/sn_final_h5ad/scvi/sn_TMCB_integration_rawcount_v2"
fn="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice/ONONH_all_rawcount"
cluster_key="majorclass"
target_cells=5000

#out="/dfs3b/ruic20_lab/junw42/HCA_TMCB/data/6_downsample/sn_TMCB_integration_rawcount_v2_downsample2k_subtype"
out="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_all_rawcount_5k"
label="majorclass"
celltype="major"

sh sb.sh -c ${comd}.sh -m 200 -t 0-1 -p fg -- $fn $cluster_key $target_cells $out $label $celltype
