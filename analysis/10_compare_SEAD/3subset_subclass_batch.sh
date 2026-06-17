#!/bin/sh

comd=HCA_ON/scripts/10_compare_SEAD/3subset_subclass

mkdir $comd

cell=Astrocyte

h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2.h5ad

sbatch -p standard --mem=200GB --time=0-1 --error=${comd}/${cell}.err --output=${comd}/${cell}.out  HCA_ON/scripts/10_compare_SEAD/3subset_subclass.sh $cell $h5ad


cell=Microglia

h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Microglia/clean/Microglia_subclass_sb.h5ad

sbatch -p standard --mem=200GB --time=0-1 --error=${comd}/${cell}.err --output=${comd}/${cell}.out  HCA_ON/scripts/10_compare_SEAD/3subset_subclass.sh $cell $h5ad


cell=Oligodendrocyte

h5ad=/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new.h5ad

sbatch -p standard --mem=200GB --time=0-1 --error=${comd}/${cell}.err --output=${comd}/${cell}.out  HCA_ON/scripts/10_compare_SEAD/3subset_subclass.sh $cell $h5ad



