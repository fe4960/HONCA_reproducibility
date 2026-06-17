#!/bin/sh
#comd="DRG/scripts/4_cellchat/1cellchat1_general"
#comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/12_CellCellInteraction/2cellchat/2cellchat_compare_ON_plot
comd=HCA_ON/scripts/12_CellCellInteraction/2cellchat/4plot_astro
#comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/12_CellCellInteraction/2cellchat/2cellchat_compare
mkdir $comd

#wd="/dfs3b/ruic20_lab/junw42/DRG/data/4_cellchat/g3_sham_downsample/"
#wd1="/dfs3b/ruic20_lab/junw42/DRG/data/4_cellchat/g3_chip_downsample/"
wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_12/"
wd1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_20to40/"
g="age_20to40_12"


sh sb.sh -c ${comd}.sh -m 60 -p s -n hpc3-l18-00 -t 0-3 -e $comd -j ${g} -- $wd $wd1 $g


wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_12/"
wd1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_70/"
g="age_70_12"


sh sb.sh -c ${comd}.sh -m 60 -p s -n hpc3-l18-00 -t 0-3 -e $comd -j ${g} -- $wd $wd1 $g



wd="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_20to40/"
wd1="/dfs3b/ruic20_lab/junw42/HCA_ON/data/11_cellchat/ONH_age_70/"
g="age_70_20to40"


sh sb.sh -c ${comd}.sh -m 60 -p s -n hpc3-l18-00 -t 0-3 -e $comd -j ${g} -- $wd $wd1 $g

