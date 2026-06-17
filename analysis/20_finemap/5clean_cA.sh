#!/bin/sh

awk '{if(($1!="NA")&&($2!="NA")){print}}' /dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq > /dfs3b/ruic20_lab/junw42/HCA_ON/data/12_snATAC/majorclass/g2p_cA_final_major_full_max_LCRE_uniq_clean
