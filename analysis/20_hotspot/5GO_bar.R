library(ggplot2)

go="enrichr_MSigDB_Hallmark_2020"

#dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_permut_hs_min_gene_200_hvg_5000_age_m"
dir1="HCA_ON/data/20_hotspot/OLIGO2_LRRC7+_permute_metacell_200_hs_min_gene_350_hvg_10000_age_m"

for( i in seq(1,7,1)){
    file=paste0(dir1,"_",i,"_bg/",go,".txt")
    df=read.table(file,header=T,sep="\t")
    df$Log10P=-log(df$P.value,base=10)
    df=df[c(1:5),]
    df=df[order(df$Log10P,decreasing=F),]
    p=ggplot(data=df,aes(x=factor(Term,levels=df$Term), y=Log10P, fill=as.numeric(Log10P)))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P-value",x="GO Term") + theme_minimal(base_size=20)+scale_fill_gradient(low = "lightblue", high = "darkblue")
#    p=ggplot(data=df,aes(x=factor(Term,levels=df$Term), y=Log10P, fill=as.numeric(Log10P)))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P-value",x="GO Term") + theme_minimal(base_size=20)+scale_fill_gradient(low = "thistle", high = "purple4")
#    out=paste0(dir1,"_",i,"_bg/",go,"_",i,"_purple.pdf")

    out=paste0(dir1,"_",i,"_bg/",go,"_",i,"_blue.pdf")
    pdf(out,height=3,width=10)
    print(p)
    dev.off()
}
