#!/bin/sh

comd=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/1gene_count_DEG.sh

err=/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/7_DEG/1gene_count_DEG

mkdir $err

#cell=Astrocyte
label=celltype
main=/dfs3b/ruic20_lab/junw42/HCA_ON/data/
dir=${main}/7_DEG/genexp_sample_subclass
mkdir -p ${dir}
ct_list="none"
for cell in Oligodendrocyte Fibroblast Endothelial_cell Astrocyte 
do
dir1=${dir}/${cell}

path=${main}/5_refine_major/scvi/lattice/ONONH_${cell}_rawcount.h5ad
list1=${main}/7_DEG/genexp_sample_subclass/snRNA_sample_num_cell_${cell}
list2=${main}/7_DEG/genexp_sample_subclass/snRNA_sample_num_${cell}

sbatch --mem=10GB --time=0-10 --error=${err}/${cell}.err --output=${err}/${cell}.out $comd ${dir1} $path ${list1} ${list2} $label ${ct_list}
done
cell=all
ct_list=/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/ct_list_major
label=majorclass
dir1=${dir}/${cell}

path=${main}/5_refine_major/scvi/lattice/ONONH_${cell}_rawcount.h5ad
list1=${main}/7_DEG/genexp_sample_subclass/snRNA_sample_num_cell_${cell}
list2=${main}/7_DEG/genexp_sample_subclass/snRNA_sample_num_${cell}

sbatch --mem=30GB --time=0-10 --error=${err}/${cell}.err --output=${err}/${cell}.out $comd ${dir1} $path ${list1} ${list2} $label ${ct_list}
