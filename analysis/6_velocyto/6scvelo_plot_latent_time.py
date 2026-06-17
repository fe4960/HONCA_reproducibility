import scanpy as sc
from matplotlib import rcParams
rcParams["figure.figsize"] = (5,5)
indir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/Astrocyte/clean/"
h5ad="Astrocyte_subclass_new_update"

hvg=2000
bname=f"{h5ad}_hvg_{hvg}_raw"

adata_veloc=sc.read(f'{indir}/veloc/{bname}_velocity.h5ad')


sc.pl.violin(adata_veloc,"latent_time",groupby="subclass")

h5ad1="Astrocyte_subclass_new2"
adata=sc.read(f'{indir}/{h5ad1}.h5ad')
adata_veloc.obs["subclass1"]=adata.obs.loc[adata_veloc.obs.index,"subclass"]


sc.settings.figdir=f"{indir}/veloc/"
sc.pl.violin(adata_veloc,"latent_time",groupby="subclass1",palette="tab20",order=["Astro_ON","Astro_ON_AFF3+","Astro_ONONH","Astro_ONH","Astro_retina","Astro_retina_NLGN1+"],rotation=90,save="_latent_time_subclass_new2.png")

sc.pl.violin(adata_veloc,"velocity_pseudotime",groupby="subclass1",palette="tab20",order=["Astro_ON","Astro_ON_AFF3+","Astro_ONONH","Astro_ONH","Astro_retina","Astro_retina_NLGN1+"],rotation=90,save="_pseudo_time_subclass_new2.png")
