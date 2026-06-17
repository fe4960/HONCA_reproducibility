library(ggplot2)
library(dplyr)
eng="enrichr_MSigDB_Hallmark_2020"

#dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_permut_hs_min_gene_200_hvg_5000_age_m"
dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_fullcell_10kgene_hs_min_gene_200_hvg_10000_scvi_m"
#dir1="HCA_ON/data/20_hotspot/Astrocyte_hs_min_gene_100_hvg_5000_scvi_m"
#dir1="HCA_ON/data/20_hotspot/Oligodendrocyte_fullcell_5kgene_hs_min_gene_250_hvg_5000_scvi_m"

df_all=NULL

go=NULL

for( i in seq(1,7,1)){
    file=paste0(dir1,"_",i,"_bg/",eng,".txt")
    df=read.table(file,header=T,sep="\t")
#    df$Log10P=-log(df$Adjusted.P.value, base=10)

        df$Log10P=-log(df$P.value,base=10)
    df1=df[c(1:5),]
#    go=rbind(go, df1$Term)    
    df1=df1[order(df1$Log10P,decreasing=F),]

    df1$module=paste0("Module ",i)
    df1$sig="darkgreen"
    df1[df1$Log10P>=2,]$sig="yellow"
#    df1=df[df$Log10P>=2,]
#    df_all=rbind(df_all,df)
#}


#df_all=df_all[(df_all$Term %in% go),]
p=ggplot(data=df1,aes(x=factor(Term,levels=df1$Term), y=Log10P, fill=sig))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P-value",x="GO Term") + theme_minimal(base_size=20)+ggtitle(unique(df1$module))#scale_fill_viridis_c(option = "viridis")+ggtitle(unique(df1$module))

#p=ggplot(data=df1,aes(x=factor(Term,levels=df$Term), y=Log10P, fill=as.numeric(Log10P)))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P.adj",x="GO Term") + theme_minimal(base_size=20)+scale_fill_viridis_c(option = "viridis")+ggtitle(unique(df$module))
#go=unique(df$module)
#out=paste0(dir1,"_",i,"_bg/module_",i,"_viridis_P_adj.pdf")

out=paste0(dir1,"_",i,"_bg/module_",i,"_viridis_reform.pdf")
n=length(df1$module)/5*3
pdf(out,height=n,width=10)
print(p)
dev.off()

}

#p <- ggplot(df1, aes(x=module,y=Term)) + geom_point(aes(size=Log10P,color=Log10P))+theme(
#    axis.line=element_line(colour="black")
#    ,panel.background = element_rect(fill="transparent") 
#    , plot.background = element_blank()
#    , panel.grid.major = element_blank()
#    , panel.grid.minor = element_blank()
#    , text=element_text(size=20))+xlab("")+ylab("")+scale_size_continuous(name=expression(paste('-',log[10],'P')),range=c(1,10),limits = c(0, 10))+scale_color_viridis_c(option = "viridis", limits = c(0, 10))


#+scale_color_gradient2( low = "blue",   mid = "white",   high = "red",   midpoint = 5, limits = c(0, 10))


#+scale_color_viridis_c(option = "viridis", limits = c(0, 10))

#+scale_fill_gradient2(  low = "blue",   mid = "white",   high = "red",   midpoint = 5, limits = c(0, 10))

#scale_fill_viridis_c(option = "viridis", limits = c(0, 10))

  
#pdf("HCA_ON/data/20_hotspot/Oligodendrocyte_precursor_cell_fullcell_10kgene_hs_min_gene_200_hvg_10000_scvi_m_dotplot.pdf",width=12)
#print(p)
#dev.off()


#    p=ggplot(data=df_all,aes(x=factor(Term,levels=df$Term), y=Log10P, fill=as.numeric(Log10P)))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P-value",x="GO Term") + theme_minimal(base_size=20)+scale_fill_viridis_c(option = "viridis")+facet_wrap(~module,ncol=4)
#    p=ggplot(data=df,aes(x=factor(Term,levels=df$Term), y=Log10P, fill=as.numeric(Log10P)))+ geom_bar(stat = "identity") + coord_flip() + labs(y="-Log10P-value",x="GO Term") + theme_minimal(base_size=20)+scale_fill_viridis_c(option = "viridis")+ggtitle(unique(df$module))
#    out=paste0(dir1,"_",i,"_bg/",go,"_",i,"_viridis.pdf")
#    out=paste0(dir1,"_all_bg_",go,"_all_blue.pdf")
#    pdf(out,height=3,width=10)
#    print(p)
#    dev.off()
