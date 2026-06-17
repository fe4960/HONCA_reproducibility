import pandas as pd
dir1="/dfs3b/ruic20_lab/junw42/gRNA_dropoff/data/sujung/scRNA/250510_X/cellranger/crispr_analysis"
singlet=pd.read_csv(f"{dir1}/singlet_list",header=0)

import scanpy as sc

dir="/dfs3b/ruic20_lab/junw42"
adata=sc.read(f"{dir}/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte_precursor_cell/clean/Oligodendrocyte_precursor_cell_hvg10000_epochnone_seurat_v3_rs_0.3_clean_res0.3rm1_sb_none_scvi_trg.h5ad")

mk={"other":["GAP43","LAMA2","EBF1","FOXP1"],"OPC":["PDGFRA","CSPG4","OLIG2"],"COP":["GPR17","VCAN","NEU4","BMP4","FYN","SOX6","PTPRZ1","PDCD4"],"NFOL":["TCF7L2","CEMIP2"],"MFOL":["MBP","OPALIN","CTPS1"],"MOL":["MBP","PLP1","MOG","GRM3","KLK6","PTGDS","ANXA5","HOPX"]}


indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/veloc"
bname="Oligodendrocyte_subclass_seurat_hvg_2000_raw"
adata=sc.read(f'{indir}/{bname}_velocity.h5ad')

celltype="Oligodendrocyte"
indir="/dfs3b/ruic20_lab/junw42"
outdir=f"{indir}/HCA_ON/data/5_refine_major/scvi/{celltype}/clean/"
sc.settings.figdir = f"{outdir}/figures"
rcParams["figure.figsize"] = (5,5)
#sc.pl.violin(adata,"latent_time",groupby="subclass", rotation=90, order=["OPC_GLIS3+","OPC","COP_NFOL","COP_LAMA2+","MFOL"],save="pseudo_time_violin.pdf")
sc.pl.violin(adata,"latent_time",groupby="subclass", order=["OLIGO_LRRC7+","OLIGO_SVEP1+","OLIGO_non_specific","OLIGO_RBFOX1+"], rotation=90,save="pseudo_time_violin.pdf")

