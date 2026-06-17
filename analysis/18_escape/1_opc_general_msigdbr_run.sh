#!/bin/sh
err=HCA_ON/scripts/18_escape/1_opc_general_msigdbr
mkdir $err

#name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat"
name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_subclass_5k_simple"
c="Oligodendrocyte_5k"
sh sb.sh -c ${err}.sh -e $err -j $c -m 50 -t 0-22 -p fg  -- $name $c

#/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_subclass_5k_simple.rds
name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_subclass_5k_simple"
c="OPC_5k"
sh sb.sh -c ${err}.sh -e $err -j $c -m 50 -t 0-22 -p fg  -- $name $c


name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new4_clean"
c="Astrocyte"
#sh sb.sh -c HCA_ON/scripts/18_escape/1_opc_general.sh -e $err -j $c -m 50 -t 0-5 -p fg -- $name $c

name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new"
c="Fibroblast"
#sh sb.sh -c HCA_ON/scripts/18_escape/1_opc_general.sh -e $err -j $c -m 50 -t 0-5 -p fg -- $name $c


name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_simple"
c="Astrocyte"
sh sb.sh -c ${err}.sh -e $err -j $c -m 50 -t 0-22 -p fg -- $name $c


name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_5k_simple"
c="Astrocyte_5k"
sh sb.sh -c ${err}.sh -e $err -j $c -m 50 -t 0-22 -p fg -- $name $c


name="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location_subclass_5k_simple"
c="Fibroblast_5k"
sh sb.sh -c ${err}.sh -e $err -j $c -m 50 -t 0-22 -p fg -- $name $c

