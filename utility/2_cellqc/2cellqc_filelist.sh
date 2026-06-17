#!/bin/sh
#ls /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE267301/cellranger/*_rename/*_rename/outs | grep 1_cellranger | awk '{split($1,a,"/"); split(a[11],b,"_"); print b[1]"\t"$_}' | sed -e "s/\://g" > /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE267301/cellranger/cellranger_list_full


#ls /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/*/*/outs | grep 1_cellranger | awk '{split($1,a,"/"); split(a[11],b,"_"); print b[1]"\t"$_}' | sed -e "s/\://g" > /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/cellranger_list_full


ls /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/*/*/outs | grep 1_cellranger | awk '{split($1,a,"/"); split(a[11],b,"_"); print b[1]"\t"$_}' | sed -e "s/\://g" > /dfs3b/ruic20_lab/junw42/HCA_ON/data/1_cellranger/pub/GSE167494/cellranger_sam/cellranger_list_full1
