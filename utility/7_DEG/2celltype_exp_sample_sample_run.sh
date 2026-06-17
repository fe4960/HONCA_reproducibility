#!/bin/sh
comd=HCA_ON/scripts/7_DEG/2celltype_exp_sample_sample
mkdir ${comd}
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/" #genexp_sample_raw_atlas_final"
od="genexp_sample_raw_atlas_final"
#main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass"
#for cell in all  #Astrocyte  Endothelial_cell  Oligodendrocyte Fibroblast 
for cell in major #all #Endothelial_cell
do	
sample_list=${main}/snRNA_sample_majorclass_num_atlas_final

#sample_list=${main}/snRNA_sample_num_${cell}
ct_list=${main}/ct_list_${cell}
indir=${main}/${od}

#indir=${main}/${cell}
#outdir=${main}/${cell}_subclass
outdir=${main}/${cell}_majorclass_sample

sbatch --mem=10GB --time=0-2 --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh ${sample_list} ${ct_list} $indir $outdir
done
