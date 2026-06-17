#!/bin/sh
ls HCA_ON/data/2_cellqc/pub/GSE267301/result/*h5ad | awk -F "/" '{print $NF}' | sed -e "s/\.h5ad//g"  | awk '{print $1"\t/dfs3b/ruic20_lab/junw42/HCA_ON/data/2_cellqc/pub/GSE267301/result/"$1".h5ad"}' > HCA_ON/data/2_cellqc/pub/GSE267301_h5ad

ls HCA_ON/data/2_cellqc/pub/GSE167494/result/*h5ad | awk -F "/" '{print $NF}' | sed -e "s/\.h5ad//g"  | awk '{print $1"\t/dfs3b/ruic20_lab/junw42/HCA_ON/data/2_cellqc/pub/GSE167494/result/"$1".h5ad"}' > HCA_ON/data/2_cellqc/pub/GSE167494_h5ad
