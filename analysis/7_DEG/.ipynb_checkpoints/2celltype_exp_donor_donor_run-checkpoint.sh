#!/bin/sh
comd=HCA_ON/scripts/7_DEG/2celltype_exp_donor_donor
mkdir ${comd}
main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass_final_update"

#main="/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_donor_subclass"
#for cell in all  #Astrocyte  Endothelial_cell  Oligodendrocyte Fibroblast 
for cell in Astrocyte Oligodendrocyte Oligodendrocyte_precursor_cell Endothelial_cell Fibroblast Microglia RGC MG  Macrophage #all #Endothelial_cell
do	
sample_list=${main}/snRNA_donor_num_${cell}
ct_list=${main}/ct_list_${cell}
indir=${main}/${cell}
outdir=${main}/${cell}_subclass
sbatch -p free --mem=10GB --time=0-2 --error=${comd}/${cell}.err --output=${comd}/${cell}.out ${comd}.sh ${sample_list} ${ct_list} $indir $outdir
done
