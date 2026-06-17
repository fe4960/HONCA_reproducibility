#!/bin/sh

dir=/dfs3b/ruic20_lab/junw42/software/deepsas/

cd $dir
source .venv/bin/activate
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_seurat_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/deepsas"
exp="Oligodendrocyte_subclass_seurat"
mkdir $outdir


#uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples


h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new5_clean_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/deepsas"
exp="Astrocyte_subclass_new5_clean"
mkdir $outdir


#uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples


h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/deepsas"
exp="Fibroblast_subclass_clean_new_rename_rmRPE_seurat_v3_rename_location"
mkdir $outdir


#uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples


h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_subclass_seurat_cycling_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/deepsas"
exp="Oligodendrocyte_precursor_cell_subclass_seurat_cycling"
mkdir $outdir

ct="Endothelial_cell"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_rmRPE_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas"
exp="${ct}_subclass_rmRPE"
mkdir $outdir


uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples



ct="Mural_cell"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_rmRPE_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas"
exp="${ct}_subclass_rmRPE"
mkdir $outdir


uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples


ct="Microglia"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_sb_clean_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas"
exp="${ct}_subclass_rmRPE"
mkdir $outdir


uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples



ct="Macrophage"
h5ad="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/${ct}_subclass_sb_clean_rmRPE_clusters.h5ad"
outdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/${ct}/clean/deepsas"
exp="${ct}_subclass_rmRPE"
mkdir $outdir


uv run python subsampling.py --input $h5ad --subsample_size 30000 --output ${outdir}/subsamples

