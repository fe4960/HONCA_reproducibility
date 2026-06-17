#!/bin/sh
err="/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/3_scanpy_check/0scanpy_check_batch"
mkdir $err
celltype="major"
in_dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellqc/pri/result/"

for queryname in BCM_22_0101_ON_RNA  BCM_22_0524_ON_RNA  BCM_22_0849_ON_RNA   BCM_22_0896_ON_RNA BCM_22_0200_ON_RNA  BCM_22_0784_ON_RNA  BCM_22_0890_ONH_RNA  BCM_23_0358_ON_RNA 
do
querylist="/dfs3b/ruic20_lab/junw42/HCA_ON/data/cellqc/pri/reaction_sample/${queryname}"
sbatch -p standard --mem=20GB --time=4-0 --account=ruic20_lab --output=${err}/${queryname}.out --error=${err}/${queryname}.err /dfs3b/ruic20_lab/junw42/HCA_ON/scripts/3_scanpy_check/0scanpy_check_batch.sh ${querylist} ${celltype} $queryname ${in_dir}
done
