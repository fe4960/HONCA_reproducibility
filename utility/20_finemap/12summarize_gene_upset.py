import sys

label_sign={}
sep="&"
#fn=["Currant_36848389PG_IST","Currant_36848389PG_ONL","Currant_36848389PG_OST","FritscheLG2016_26691988NG","Gharahkhani_33627673NC","Hysi_32231278NG","Jiang_34737426NG"]
fn=["Khawaja2018_29785010IOP", "Gharahkhani_33627673NC", "Khor2016_27064256PCAG", "Springelkamp2017_28073927CA", "Springelkamp2017_28073927DA", "Springelkamp2017_28073927IOP", "Springelkamp2017_28073927VCDR"]
for fn1 in fn:
#	file1=f'/storage/chenlab/Users/junwang/human_meta/data/finemap/{fn1}_pip_var_all_merged_uniform_anno'
    file1=f'/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/{fn1}_pip_var_all_merged_uniform_anno_ON'

#	file1=f'/storage/chenlab/Users/junwang/human_meta/data/age_DESeq2_batch_age_dream/{ct}_interval_cpm01_snRNA_rmH/{ct}_DEG_res'
    with open (file1,"r") as f1:
        for line in f1:
            info=line.rstrip().split()
            if (info[-4] != "anno;") and ( float(info[-5]) >=0.1):
				#anno;3_UTR|peak|
#				feature=info[-2].replace("anno;","").replace("3_UTR","Exon").replace("5_UTR","Exon").split("|")
                feature=info[-4].replace("anno;","").split("|")
                ft_line=sep.join(feature)
                if ft_line not in label_sign:
                    label_sign[ft_line]={}
                label_sign[ft_line][info[1]]=1



#output=f'/storage/chenlab/Users/junwang/human_meta/data/finemap/pip_var_all_merged_uniform_anno_all'
#output=f'/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/pip_var_all_merged_uniform_anno_all_ON_02'
output=f'/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/pip_var_all_merged_uniform_anno_all_ON_01'

#output=f'/dfs3b/ruic20_lab/junw42/reference/GWAS/finemap_sumstat/suise/pip_var_all_merged_uniform_anno_all_ON_05'
#output=f'/storage/chenlab/Users/junwang/human_meta/data/finemap/pip_var_all_merged_uniform_anno_all_pip01'

out=open(output,"w") 

arr=[]
for label in label_sign:
    num=len(label_sign[label].keys())
   arr.append(f'"{label}"={num}')

#	out.write(f'"{label}"={num},')
#    arr.append(f'"{label}"={num}')

new_line="\n".join(arr)

out.write(f'{new_line}')

out.close()

