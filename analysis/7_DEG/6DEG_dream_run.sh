#!/bin/sh

comd=HCA_ON/scripts/7_DEG/6DEG_dream

#awk -F "," '{a=substr($1,1,3); print $_","a}' /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race > /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch
#cut -d "," -f 2-7,9- /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch | sort -u > /dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid

mkdir $comd

#for c in  all #Oligodendrocyte   # 
#for c in  all #Fibroblast   #Astrocyte  #Endothelial_cell  #all Astrocyte  Endothelial_cell  Fibroblast
od="genexp_donor_subclass_final"
do
#file=HCA_ON/data/7_DEG/genexp_sample_subclass/ct_list_${c}
file=HCA_ON/data/7_DEG/${od}/ct_list_${c}

while read line
do
ct=$line
sbatch --nodelist=hpc3-l18-05 --account=ruic20_lab -p standard --mem=20GB --time=0-2 --error=${comd}/${line}.err --output=${comd}/${line}.out ${comd}.sh $c $ct
done < $file
done

