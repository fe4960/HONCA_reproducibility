import scanpy as sc
import anndata as ad
rpe=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/RPEChoroid/Lattice/RPEChoroid_h5ad/snRNA_fixsampleid/RPEChoroid_snRNA_allcells_rawcounts.h5ad")
fb_rpe_ad=rpe[rpe.obs["majorclass"]=="RPE_Adult"]
fb_rpe_ad.obs["celltype"]=fb_rpe_ad.obs["subclass"]
fb_rpe_ft=rpe[rpe.obs["majorclass"]=="RPE_Fetal"]
fb_rpe_ft.obs["celltype"]=fb_rpe_ft.obs["subclass"]

fb_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/RPE_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#fb_on=fb_on[fb_on.obs["leiden"]!='6']
fb_on.X=fb_on.layers["counts"]
fb_on.obs["celltype"]=fb_on.obs["leiden"]
adata_list=[fb_rpe_ad,fb_rpe_ft , fb_on]
adata=ad.concat(adata_list)
adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/RPE/clean/RPE_rpe_osf_tmcb_ononh.h5ad")
