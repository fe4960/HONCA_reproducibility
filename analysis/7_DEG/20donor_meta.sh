#!/bin/sh
awk -F "," '{print $2"\t"$3"\t"$4"\t"$7}' /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid  | sort -u > /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor
