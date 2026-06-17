import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#trait=["raceAfrican", "raceEur", "raceAdmixed", "age", "gender", "ageGender","tissueFovea","tissueMacular"]
#trait=["tissueFovea","tissueMacular"]
trait=["age"]
annot_kws={'fontsize': 11 }
#cell=["Oligodendrocyte"]
cell=["all"]
#cell=["Endothelial_cell","Astrocyte","Fibroblast"]
for tr in trait:
    for ct in cell:
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_eff_ct"
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_eff_pair"
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_eff"
        fn=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{tr}_DEG_res_cpm_{tr}_cell_eff_ct_ancFix"
        cor=pd.read_csv(fn,sep="\t",header=0)

        matrix = np.triu(np.ones_like(cor))

        df = pd.DataFrame(matrix)
        df.columns=cor.columns
        df.index=cor.columns
        cor.index=cor.columns
        cor=np.multiply(cor,df)
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_dir_ct"
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_dir_pair"
#	fn=f"/storage/chenlab/Users/junwang/human_ret_anc/data/DEG/region_race_batch_ageNum_dream_atlas_sample/{tr}_DEG_res_cpm_{tr}_cell_dir"
        fn=f"/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/{ct}_subclass_DEG/{tr}_DEG_res_cpm_{tr}_cell_dir_ct_ancFix"
        cor1=pd.read_csv(fn,sep="\t",header=0)

        matrix1 = np.tril(np.ones_like(cor))
        df1 = pd.DataFrame(matrix1)

        df1.columns=cor1.columns
        df1.index=cor1.columns
        cor1.index=cor1.columns
        cor1=np.multiply(cor1,df1)

        cor2=cor+cor1

        np.fill_diagonal(cor2.values, 1)

	# Add row names
#	df.index = cor.columns

        plt.figure(figsize=(12,12))
        plt.title(f"{tr} effect size")
	# using the upper triangle matrix as mask
        g=sns.heatmap(cor2, annot=True,  annot_kws=annot_kws )
        g.set_xticklabels(g.get_xticklabels(), rotation=90)
        plt.savefig(f'{fn}_eff.pdf',dpi=500)
        plt.close("all") 

