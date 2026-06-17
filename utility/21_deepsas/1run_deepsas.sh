#!/bin/sh



dir=/dfs3b/ruic20_lab/junw42/software/deepsas/

cd $dir
source .venv/bin/activate
#h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_clusters.h5ad"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/deepsas"
#exp="Oligodendrocyte_subclass_seurat"
exp="Data_0"
mkdir $outdir
uv run python -u ${dir}/deepsas_v1.py --seed 7  --input_data_count $h5ad  --output_dir $outdir   --exp_name exp --device_index 0 --retrain #> ./example.log
#uv run python -u ${dir}/deepsas_v1.py --seed 7 --batch_id sampleid --input_data_count $h5ad  --output_dir $outdir   --exp_name exp --device_index 0 --retrain > ./example.log
