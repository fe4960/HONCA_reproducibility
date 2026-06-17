import scanpy as sc
import pandas as pd

dir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/lattice"


meta=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/scripts/5_refine_major/Table_HCA_donorinfo.txt_format_uniq",index_col=0,header=0,sep="\t")

pmit=meta.to_dict()["pmit"]

globe=meta.to_dict()["globe"]

meta_sanes=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/3_ref/joshSanes_meta/donor_info_pmit_odos",index_col=0,header=0,sep="\t")

#pmits=meta_sanes.to_dict()["pmit"]

#globes=meta_sanes.to_dict()["globe"]



def correct_sample(adata):

    s1="BCM_22_0890_ONH_RNA_s2"
    s2="BCM_22_0458_ONH_RNA_s2"
    adata.obs["sampleid"]=adata.obs["sampleid"].cat.add_categories([s2])
    adata.obs["reactionid"]=adata.obs["reactionid"].cat.add_categories([s2])
    t="25"
    if t not in adata.obs["age"].values:
        adata.obs["age"]=adata.obs["age"].cat.add_categories(["25"])

    t="BCM_22_0458"
    if t not in adata.obs["donor"].values:
        adata.obs["donor"]=adata.obs["donor"].cat.add_categories(["BCM_22_0458"])
#    if 25 not in adata.obs["age_year"].values:
#        adata.obs["age_year"]=adata.obs["age_year"].cat.add_categories([25])



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
    t="MMD_23_20181"
    if t in adata.obs["donor"].values:
        print("in")
    else:
        adata.obs["donor"]=adata.obs["donor"].cat.add_categories(["MMD_23_20181"])

    t="16"
    if t not in adata.obs["age"].values:
        adata.obs["age"]=adata.obs["age"].cat.add_categories(["16"])
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

    return adata




        
        
        
adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Fibroblast/clean/Fibroblast_subclass_LC_other1.h5ad")
adata.obs["majorclass"]=adata.obs["majorclass1"]

tag=['connectivities', 'distances']
for t in tag:
    del adata.obsp[t]


tag=["majorclass1", "leiden"]
for t in tag:
    del adata.obs[t]

tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts']
#tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm']
for t in tag:
    del adata.var[t]

tag=['hvg', 'leiden', 'leiden_colors',"rank_genes_groups",'_scvi_manager_uuid', '_scvi_uuid', 'dendrogram_subclass',  'neighbors']
for t in tag:
    del adata.uns[t]

del adata.raw
adata=correct_sample(adata)


adata.obs["pmit"]=adata.obs["donor"].copy()
adata.obs["pmit"]=adata.obs["pmit"].replace(pmit)

adata.obs["globe"]=adata.obs["donor"].copy()
adata.obs["globe"]=adata.obs["globe"].replace(globe)


for i in meta_sanes.index:
    adata.obs.loc[adata.obs["sampleid"]==i,"pmit"]=meta_sanes.loc[i,"pmit"]
    adata.obs.loc[adata.obs["sampleid"]==i,"globe"]=meta_sanes.loc[i,"globe"]



tmp=["sampleid", "donor", "reactionid", "sampleid_legacy","subclass","majorclass","globe","pmit"]
    
for t in tmp:
    adata.obs[t]=adata.obs[t].cat.remove_unused_categories()

adata.obs["celltype"]=adata.obs["subclass"].copy()
    
del adata.obs["subclass"]



adata1=adata.copy()
del adata1.layers["counts"]
adata1.write(f"{dir}/ONONH_Fibroblast_normcount.h5ad")

adata.X=adata.layers["counts"].copy()
del adata.layers["counts"]
adata.write(f"{dir}/ONONH_Fibroblast_rawcount.h5ad")



adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/Astrocyte_subclass_new2.h5ad")
adata.obs["majorclass"]=adata.obs["majorclass1"]

tag=['connectivities', 'distances']
for t in tag:
    del adata.obsp[t]

tag=["majorclass1", "leiden"]
for t in tag:
    del adata.obs[t]


#tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm']
tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts']

for t in tag:
    del adata.var[t]

tag=['hvg', 'leiden', 'leiden_colors', "rank_genes_groups",'_scvi_manager_uuid', '_scvi_uuid', 'dendrogram_subclass',   'neighbors']
for t in tag:
    del adata.uns[t]

del adata.raw
adata=correct_sample(adata)


adata.obs["pmit"]=adata.obs["donor"].copy()
adata.obs["pmit"]=adata.obs["pmit"].replace(pmit)

adata.obs["globe"]=adata.obs["donor"].copy()
adata.obs["globe"]=adata.obs["globe"].replace(globe)

for i in meta_sanes.index:
    adata.obs.loc[adata.obs["sampleid"]==i,"pmit"]=meta_sanes.loc[i,"pmit"]
    adata.obs.loc[adata.obs["sampleid"]==i,"globe"]=meta_sanes.loc[i,"globe"]


tmp=["sampleid", "donor", "reactionid", "sampleid_legacy","subclass","majorclass","globe","pmit"]

for t in tmp:
    adata.obs[t]=adata.obs[t].cat.remove_unused_categories()

adata.obs["celltype"]=adata.obs["subclass"].copy()

del adata.obs["subclass"]


adata1=adata.copy()
del adata1.layers["counts"]
adata1.write(f"{dir}/ONONH_Astrocyte_normcount.h5ad")

adata.X=adata.layers["counts"].copy()
del adata.layers["counts"]
adata.write(f"{dir}/ONONH_Astrocyte_rawcount.h5ad")



adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Endothelial_cell/clean/Endothelial_cell_subclass_sb.h5ad")
adata.obs["majorclass"]=adata.obs["majorclass1"]

tag=['connectivities', 'distances']
for t in tag:
    del adata.obsp[t]


tag=["majorclass1", "leiden"]
for t in tag:
    del adata.obs[t]

tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts']

#tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm']
for t in tag:
    del adata.var[t]

tag=['hvg', 'leiden', 'leiden_colors',"rank_genes_groups",'_scvi_manager_uuid', '_scvi_uuid', 'dendrogram_subclass',  'neighbors']
for t in tag:
    del adata.uns[t]

del adata.raw
adata=correct_sample(adata)


adata.obs["pmit"]=adata.obs["donor"].copy()
adata.obs["pmit"]=adata.obs["pmit"].replace(pmit)

adata.obs["globe"]=adata.obs["donor"].copy()
adata.obs["globe"]=adata.obs["globe"].replace(globe)

for i in meta_sanes.index:
    adata.obs.loc[adata.obs["sampleid"]==i,"pmit"]=meta_sanes.loc[i,"pmit"]
    adata.obs.loc[adata.obs["sampleid"]==i,"globe"]=meta_sanes.loc[i,"globe"]


tmp=["sampleid", "donor", "reactionid", "sampleid_legacy","subclass","majorclass","globe","pmit"]

for t in tmp:
    adata.obs[t]=adata.obs[t].cat.remove_unused_categories()

adata.obs["celltype"]=adata.obs["subclass"].copy()

del adata.obs["subclass"]


adata1=adata.copy()
del adata1.layers["counts"]
adata1.write(f"{dir}/ONONH_Endothelial_cell_normcount.h5ad")

adata.X=adata.layers["counts"].copy()
del adata.layers["counts"]
adata.write(f"{dir}/ONONH_Endothelial_cell_rawcount.h5ad")


adata=sc.read("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Oligodendrocyte/clean/Oligodendrocyte_subclass_new.h5ad")
adata.obs["majorclass"]=adata.obs["majorclass1"]


tag=['connectivities', 'distances']
for t in tag:
    del adata.obsp[t]


tag=["majorclass1", "leiden"]
for t in tag:
    del adata.obs[t]

tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts']

#tag=['highly_variable', 'means', 'dispersions', 'dispersions_norm', 'highly_variable_nbatches', 'highly_variable_intersection', 'highly_variable_rank', 'variances', 'variances_norm']
for t in tag:
    del adata.var[t]

tag=['hvg', 'leiden', 'leiden_colors', "rank_genes_groups",'_scvi_manager_uuid', '_scvi_uuid', 'dendrogram_subclass',   'neighbors']
for t in tag:
    del adata.uns[t]

del adata.raw
adata=correct_sample(adata)


adata.obs["pmit"]=adata.obs["donor"].copy()
adata.obs["pmit"]=adata.obs["pmit"].replace(pmit)

adata.obs["globe"]=adata.obs["donor"].copy()
adata.obs["globe"]=adata.obs["globe"].replace(globe)

for i in meta_sanes.index:
    adata.obs.loc[adata.obs["sampleid"]==i,"pmit"]=meta_sanes.loc[i,"pmit"]
    adata.obs.loc[adata.obs["sampleid"]==i,"globe"]=meta_sanes.loc[i,"globe"]

tmp=["sampleid", "donor", "reactionid", "sampleid_legacy","subclass","majorclass","globe","pmit"]

for t in tmp:
    adata.obs[t]=adata.obs[t].cat.remove_unused_categories()

adata.obs["celltype"]=adata.obs["subclass"].copy()

del adata.obs["subclass"]

adata1=adata.copy()
del adata1.layers["counts"]
adata1.write(f"{dir}/ONONH_Oligodendrocyte_normcount.h5ad")

adata.X=adata.layers["counts"].copy()
del adata.layers["counts"]
adata.write(f"{dir}/ONONH_Oligodendrocyte_rawcount.h5ad")

