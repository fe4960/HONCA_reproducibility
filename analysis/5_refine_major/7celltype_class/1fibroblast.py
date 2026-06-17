import scanpy as sc
import anndata as ad
#rpe=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/RPEChoroid/Lattice/RPEChoroid_h5ad/snRNA_fixsampleid/RPEChoroid_snRNA_allcells_rawcounts.h5ad")
rpe=sc.read("/dfs3b/ruic20_lab/jianms1/Single_cell/scRNA_snRNA_1-7-2025/snRNA/Stromal/annotation/snRNA_Stromal_leiden_1_0.3_0.21_subclass_annotated.h5ad")
fb_rpe=rpe[rpe.obs["majorclass_stromal"].isin(["Fibroblasts","Fibroblasts_Fetal"])]
fb_rpe.obs["celltype"]=fb_rpe.obs["subclass_stromal"]
fb_rpe.obs["origin"]="RPE"
osf=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/OcularSurface/Lattice/OcularSurface_h5ad/snRNA_fixbarcode/OcularSurface_snRNA_allcells_rawcounts.h5ad")
fb_osf=osf[osf.obs["majorclass"]=="Fibroblasts"]
fb_osf.obs["majorclass"]="Fibroblast"
fb_osf.obs["origin"]="OcularSurface"
tmcb=sc.read("/dfs3b/ruic20_lab/jinl14/mrrdir/wkfl/hrcaprj/HECA/TMCB/Lattice/TMCB_h5ad/TMCB_snRNA_allcells_rawcounts.h5ad")
fb_tmcb=tmcb[tmcb.obs["majorclass"]=="Fibroblast"]
fb_tmcb.obs["origin"]="TMCB"
#fb_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg2000_epochnone_seurat_v3_rs_0.3_clean_rm_scvi_trg.h5ad")
#fb_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_LC_other1.h5ad")
fb_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_clean_new.h5ad")

#fb_on=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_hvg2000_epochnone_seurat_v3_rs_0.3_clean_scvi_trg.h5ad")
#fb_on=fb_on[fb_on.obs["leiden"]!='6']
fb_on.X=fb_on.layers["counts"]
fb_on.obs["celltype"]=fb_on.obs["subclass"]
fb_on.obs["origin"]="ONONH"
#fb_on.obs["celltype"]=fb_on.obs["leiden"]
adata_list=[fb_rpe, fb_osf, fb_tmcb, fb_on]
adata=ad.concat(adata_list)
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh.h5ad")
#adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh_clean.h5ad")

adata.write("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_rpe_osf_tmcb_ononh_clean_final_update.h5ad")
