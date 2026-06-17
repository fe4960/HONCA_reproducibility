#!/bin/sh
comd=HCA_ON/scripts/7_DEG/14GO_enrichR_sub

od="NA"
mkdir $comd
dem="age_cpm1"
#for cell in Oligodendrocyte 
for cell in Fibroblast #Astrocyte
do
#file=HCA_ON/data/7_DEG/ct_list_${cell}
#file=HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/ct_list_${cell}3
file=HCA_ON/data/7_DEG/genexp_donor_subclass_final_new/ct_list_${cell}

j=1
while read line
do
#j=$(($j+1))
ct=$line
sbatch  --account=ruic20_lab -p free-gpu --mem=10GB --time=0-1 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct  $od $dem

done < $file
done




dem="age_range_30_60_cpm5"

for cell in major 
do
file=HCA_ON/data/7_DEG/ct_list_${cell}

j=1
while read line
do
#j=$(($j+1))
ct=$line

#sbatch  --account=ruic20_lab -p free-gpu --mem=10GB --time=0-1 --error=${comd}/${cell}_${ct}.err --output=${comd}/${cell}_${ct}.out ${comd}.sh $cell $ct  $od $dem

done < $file
done
