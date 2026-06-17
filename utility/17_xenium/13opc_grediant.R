library(ggplot2)
data=read.table("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_6_cell_cluster_dist_y_comb_w1.txt", sep="\t",header=T)
data1=data[data$harmony_anno=="Oligodendrocyte_precursor_cell",]


library(ggplot2)
library(ggpubr)

p=ggplot(data1, aes(x = bin_y, y = prop, color=harmony_anno)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x) + stat_cor(size=5, label.y=0.015) + xlab("Bin") + ylab("Cell proportion")+ facet_wrap(~harmony_anno, scale="free_y") +
    theme_classic() +
    theme(
        legend.position="none",
        axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 25),
    axis.text.y = element_text(size = 25)
    )
#, label.y=0
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/ONONH_comb_6_cell_cluster_dist_y_comb_w1_opc.pdf", width=6.5, height=5)

#pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/RNASEH2B_DEG.pdf", width=7, height=6)
print(p)
dev.off()

