import scanpy as sc
from matplotlib import rcParams
import pandas as pd
rcParams["figure.figsize"] = (5,5)
sc.settings.figdir = "HCA_ON/data/5_refine_major/scvi/major/figures"
fl=["hvg10000_sb_epochnone_seurat_v3_rs"]

meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt_format_uniq",index_col=0,header=0,sep="\t")

pmit=meta.to_dict()["pmit"]

globe=meta.to_dict()["globe"]

meta_sanes=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/donor_info_pmit_odos",index_col=0,header=0,sep="\t")

for f1 in fl:
    adata=sc.read(f"HCA_ON/data/5_refine_major/scvi/ON_ONH_clean_new_update_final.h5ad")
    del adata.obsm["X_umap"]
    del adata.obsm["X_scVI"]


    adata.obs["pmit"]=adata.obs["donor"].copy()
    adata.obs["pmit"]=adata.obs["pmit"].replace(pmit)

    adata.obs["globe"]=adata.obs["donor"].copy()
    adata.obs["globe"]=adata.obs["globe"].replace(globe)


    for i in meta_sanes.index :
        adata.obs.loc[adata.obs["sampleid"]==i,"pmit"]=meta_sanes.loc[i,"pmit"]
        adata.obs.loc[adata.obs["sampleid"]==i,"globe"]=meta_sanes.loc[i,"globe"]


#########    adata1=sc.read(f"HCA_ON/data/5_refine_major/scvi/major/major_{f1}_1_clean_ON_ONH_new_scvi_trg.h5ad")
#######    adata=adata[adata1.obs.index]
#    adata.uns["umap"]=adata1.uns["umap"]
#    adata.uns["neighbors"]=adata1.uns["neighbors"]
#######    adata.obsm["X_umap"]=adata1.obsm["X_umap"]
#######    adata.obsm["X_scVI"]=adata1.obsm["X_scVI"]
#    adata.obsp["connectivities"]=adata1.obsp["connectivities"]
#    adata.obsp["distances"]=adata1.obsp["distances"]

#    adata=adata[~adata.obs.index.isin(idx)]
#    adata.obs["subclass"]=obs.loc[adata.obs.index,"subclass"]
    adata.obs["majorclass1"]=adata.obs["majorclass1"].cat.add_categories([ "Immune_cell", "Pigmented_cell", "Unknown?", "T_cell", "NK_cell","Mast_cell","Dendritic_cell","B_cell_activated","B_cell", "Immune_RBMS3"])
    adata.obs.loc[adata.obs["subclass"].isin(["Pigmented_cell_1","Pigmented_cell_2"]),"majorclass1"]="Pigmented_cell"
        
#    adata.obs.loc[adata.obs["majorclass1"]=="NK_T_cell","majorclass1"]=adata.obs.loc[adata.obs["majorclass1"]=="NK_T_cell","subclass"].astype(str)
    adata.obs.loc[adata.obs["majorclass1"]=="NK_T_cell","majorclass1"]="Immune_cell"

    adata.obs.loc[adata.obs["majorclass1"]=="B_cell_activated","subclass"]="B_cell"
    adata.obs.loc[adata.obs["majorclass1"]=="B_cell_activated","majorclass1"]="Immune_cell"

    adata.obs["majorclass1"]=adata.obs["majorclass1"].replace("Unknown?","Immune_RBMS3")
    
    adata.obs["majorclass1"]=adata.obs["majorclass1"].cat.remove_unused_categories()

    adata=adata[adata.obs["majorclass1"]!="Immune_RBMS3"]

#    adata.obs["subclass"]=adata.obs["subclass"].cat.add_categories([ "OLIGO_SVEP1-hi", "OLIGO_LRRC7-hi", "OLIGO_RBFOX1-hi", "OLIGO_non-specific", "OLIGO_NAV2-hi"])

#    adata.obs.loc[idx0,"subclass"]="OLIGO_RBFOX1-hi"

#    adata.obs.loc[idx1,"subclass"]="OLIGO_LRRC7-hi"

#    adata.obs.loc[idx2,"subclass"]="OLIGO_SVEP1-hi"

#    adata.obs.loc[idx3,"subclass"]="OLIGO_non-specific"

#    adata.obs.loc[idx4,"subclass"]="OLIGO_NAV2-hi"

#    adata.obs["subclass"]=adata.obs["subclass"].cat.remove_unused_categories()

    s1="BCM_22_0890_ONH_RNA_s2"
    s2="BCM_22_0458_ONH_RNA_s2"
    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="25"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=25
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Female"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="Hispanic"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_22_0458"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2

    s1="BCM_23_0229_ONH_RNA"
    s2="BCM_23_0231_ONH_RNA"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="69"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=69
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Male"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_23_0231"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2


    s1="BCM_23_0313_ON_RNA"
    s2="MMD_23_20181_ON_RNA"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="16"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=16
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Male"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="MMD_23_20181"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2


    s1="BCM_23_0358_ON_RNA_s2"
    s2="BCM_23_0231_ON_RNA_s2"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="69"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=69
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Male"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_23_0231"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2


    s1="MMD_23_17976_ONH_RNA"
    s2="BCM_23_1080_ONH_RNA_s1"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="71"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=71
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Male"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_23_1080"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2


    s1="BCM_22_0200_ONH_RNA"
    s2="BCM_22_0896_ONH_RNA_s2"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="62"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=62
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Female"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_22_0896"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2


    s1="MMD_23_20181_ONH_RNA"
    s2="BCM_23_1080_ONH_RNA_s2"

    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])

    adata.obs.loc[adata.obs["reactionid"]==s1,"age"]="71"
    adata.obs.loc[adata.obs["reactionid"]==s1,"age_year"]=71
    adata.obs.loc[adata.obs["reactionid"]==s1,"gender"]="Male"
    adata.obs.loc[adata.obs["reactionid"]==s1,"race"]="White"
    adata.obs.loc[adata.obs["reactionid"]==s1,"donor"]="BCM_23_1080"
    adata.obs.loc[adata.obs["reactionid"]==s1,"sampleid"]=s2
    adata.obs.loc[adata.obs["reactionid"]==s1,"reactionid"]=s2

    adata.obs["majorclass"]=adata.obs["majorclass1"]
############    del adata.obs["majorclass1"]
##########    del adata.obs["leiden"]

############    sc.pl.embedding(adata, basis="X_umap", color=["majorclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_ON_ONH_clean_new_majorclass.png")
###########    sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_ON_ONH_clean_new_subclass.png")


#    sc.pl.embedding(adata, basis="X_umap", color=["majorclass1"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_new_majorclass.png")
#    sc.pl.embedding(adata, basis="X_umap", color=["subclass"], ncols=1,frameon=False,palette="tab20", save=f"_major_{f1}_1_clean_ON_ONH_new_subclass.png")

    del adata.raw

    tmp=["sampleid", "donor", "reactionid", "sampleid_legacy","subclass","majorclass"]

    for t in tmp:
        adata.obs[t]=adata.obs[t].cat.remove_unused_categories()

    adata.obs["celltype"]=adata.obs["subclass"].copy()

    del adata.obs["subclass"]

    

    adata.write("HCA_ON/data/5_refine_major/scvi/major/major_ON_ONH_new_update_final.h5ad")
    df=adata.obs[['reactionid', 'sampleid', 'donor', 'race', 'gender', 'age', 'tissue','sampleid_legacy','age_year',"globe","pmit"]]
    df1=df.drop_duplicates()
    df1.to_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_update",index=False)

    dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice_update_final"

    adata.obs.to_csv(f"HCA_ON/data/5_refine_major/scvi/major/ONONH_new_update_final.obs.gz",sep="\t")


    adata=adata[-adata.obs["majorclass1"].isin(["Rod", "Cone", "AC","BC", "HC", "RGC", "MG", "RPE"])]
#    adata2=adata.copy()
#    del adata2.layers["counts"]
#    adata.write(f"{dir}/ONONH_all_normcount_only.h5ad")

#    adata.X=adata.layers["counts"]
#    del adata.layers["counts"]
    adata.write(f"{dir}/ONONH_all_raw_normcount_only.h5ad")

    adata.obs.to_csv(f"HCA_ON/data/5_refine_major/scvi/major/ONONH_new_update_final_ONONH_only.obs.gz",sep="\t")

#    adata2=adata.copy()
#    del adata2.layers["counts"]
#    adata2.write(f"{dir}/ONONH_all_normcount.h5ad")

#    adata.X=adata.layers["counts"]
#    del adata.layers["counts"]
#    adata.write(f"{dir}/ONONH_all_rawcount.h5ad")
