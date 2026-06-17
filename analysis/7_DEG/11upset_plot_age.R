library(UpSetR)
input=c('Fibro_RPEchoroid&Fibro_arachnoid&Fibro_dura&Fibro_perivascular&Fibro_pia'=11,'Fibro_RPEchoroid&Fibro_arachnoid&Fibro_dura&Fibro_perivascular&Fibro_pia&Fibro_sclera'=56,'Fibro_RPEchoroid'=361,'Fibro_RPEchoroid&Fibro_arachnoid&Fibro_dura&Fibro_perivascular&Fibro_pia&Fibro_sclera&Fibro_x'=874,'Fibro_RPEchoroid&Fibro_arachnoid&Fibro_dura&Fibro_pia&Fibro_sclera&Fibro_x'=75,'Fibro_pia'=41,'Fibro_RPEchoroid&Fibro_pia'=12,'Fibro_RPEchoroid&Fibro_dura&Fibro_perivascular&Fibro_pia&Fibro_sclera&Fibro_x'=11,'Fibro_RPEchoroid&Fibro_arachnoid&Fibro_dura&Fibro_pia&Fibro_sclera'=28)
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Fibroblast_subclass_DEG/Fibroblast_upset_list_age_all_ct.pdf",width=10)

upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "Fibroblast",
      sets.bar.color="blue"
      )
dev.off()

input=c('Artery_Arteriole&Vein_Venule&Venule_POSTN.&Venule_postcapillary'=16,'Artery_Arteriole'=22,'Artery_Arteriole&Capillary&Vein_Venule&Venule_POSTN.&Venule_postcapillary'=83,'Artery_Arteriole&Capillary&Capillary_fenestrated&Venule_postcapillary'=11,'Artery_Arteriole&Venule_POSTN.&Venule_postcapillary'=11,'Artery_Arteriole&Capillary&Vein_Venule&Venule_POSTN.'=57,'Artery_Arteriole&Capillary'=15,'Artery_Arteriole&Capillary&Capillary_fenestrated&Vein_Venule'=24,'Capillary&Capillary_fenestrated&Vein_Venule&Venule_POSTN.&Venule_postcapillary'=16,'Artery_Arteriole&Capillary&Capillary_fenestrated&Vein_Venule&Venule_POSTN.'=62,'Capillary&Vein_Venule&Venule_postcapillary'=11,'Artery_Arteriole&Capillary&Capillary_fenestrated&Venule_POSTN.&Venule_postcapillary'=16,'Vein_Venule'=19,'Capillary_fenestrated'=14,'Artery_Arteriole&Capillary&Vein_Venule&Venule_postcapillary'=38,'Capillary&Vein_Venule'=42,'Artery_Arteriole&Capillary_fenestrated&Vein_Venule&Venule_POSTN.&Venule_postcapillary'=32,'Capillary'=20,'Artery_Arteriole&Capillary&Capillary_fenestrated&Vein_Venule&Venule_postcapillary'=81,'Venule_postcapillary'=107,'Capillary&Capillary_fenestrated&Vein_Venule&Venule_postcapillary'=12,'Artery_Arteriole&Capillary&Vein_Venule'=72,'Artery_Arteriole&Capillary_fenestrated&Vein_Venule&Venule_postcapillary'=12,'Capillary&Vein_Venule&Venule_POSTN.'=22,'Artery_Arteriole&Capillary&Capillary_fenestrated&Vein_Venule&Venule_POSTN.&Venule_postcapillary'=927)
pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Endothelial_cell_subclass_DEG/Endothelial_cell_upset_list_age_all_ct.pdf",width=10)

upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "Endothelial cell",
      sets.bar.color="blue"
      )
dev.off()



input=c('Astro_ON&Astro_ONH&Astro_ONONH'=64,'Astro_ONH&Astro_ONONH&Astro_retina'=47,'Astro_ONH'=48,'Astro_retina'=188,'Astro_ON&Astro_ONH&Astro_ONONH&Astro_retina'=1078,'Astro_ON&Astro_ONONH'=42,'Astro_ON'=55,'Astro_ONH&Astro_ONONH'=13,'Astro_ON&Astro_ONONH&Astro_retina'=40,'Astro_ONONH'=21)


pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Astrocyte_subclass_DEG/Astrocyte_upset_list_age_all_ct.pdf",width=10)

upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "Astrocyte",
      sets.bar.color="blue"
      )
dev.off()


input=c('OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_non_specific'=14,'OLIGO_LRRC7.hi&OLIGO_RBFOX1.hi&OLIGO_SVEP1.hi'=22,'OLIGO_non_specific'=117,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi'=19,'OLIGO_LRRC7.hi&OLIGO_RBFOX1.hi&OLIGO_SVEP1.hi&OLIGO_non_specific'=36,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_RBFOX1.hi&OLIGO_SVEP1.hi'=41,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_SVEP1.hi&OLIGO_non_specific'=18,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_RBFOX1.hi&OLIGO_non_specific'=29,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_SVEP1.hi'=15,'OLIGO_LRRC7.hi'=86,'OLIGO_NAV2.hi'=79,'OLIGO_LRRC7.hi&OLIGO_NAV2.hi&OLIGO_RBFOX1.hi&OLIGO_SVEP1.hi&OLIGO_non_specific'=1017)

pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/Oligodendrocyte_subclass_DEG/Oligodendrocyte_upset_list_age_all_ct.pdf",width=10)

upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "Oligodendrocyte",
      sets.bar.color="blue"
      )
dev.off()

input=c('Astrocyte&Fibroblast'=18,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Mural_cell&Oligodendrocyte&Rod'=11,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&Rod'=105,'Fibroblast&Mural_cell'=15,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&RPE&Rod'=231,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Microglia&Mural_cell&Oligodendrocyte&Rod'=11,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Microglia&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&RPE&Rod'=21,'AC&BC&MG&Rod'=11,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Microglia&Mural_cell&Oligodendrocyte_precursor_cell&Rod'=28,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Oligodendrocyte&Rod'=17,'Oligodendrocyte_precursor_cell'=45,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Microglia&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&Rod'=40,'Endothelial_cell&Fibroblast&Mural_cell'=14,'Endothelial_cell&Fibroblast'=20,'Astrocyte'=42,'Endothelial_cell&Fibroblast&Oligodendrocyte'=11,'Fibroblast'=227,'Astrocyte&Endothelial_cell&Fibroblast&Mural_cell'=12,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&RPE&Rod'=12,'Fibroblast&Oligodendrocyte'=20,'Endothelial_cell'=113,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Mural_cell&Oligodendrocyte&Rod'=13,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Oligodendrocyte&RPE&Rod'=26,'Endothelial_cell&Fibroblast&Macrophage&Mural_cell&Oligodendrocyte'=10,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Mural_cell&Oligodendrocyte_precursor_cell&Rod'=11,'Oligodendrocyte'=105,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Mural_cell&Oligodendrocyte&Oligodendrocyte_precursor_cell&Rod'=11,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Oligodendrocyte_precursor_cell&Rod'=28,'Mural_cell'=240,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Macrophage&Microglia&Mural_cell&Rod'=11,'Astrocyte&Oligodendrocyte'=10,'AC&Astrocyte&BC&Endothelial_cell&Fibroblast&MG&Microglia&Mural_cell&Rod'=11,'Astrocyte&Endothelial_cell&Fibroblast'=14,'Rod'=12,'Microglia'=15)

pdf("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_subclass/all_subclass_DEG/all_upset_list_age_all_ct.pdf",width=10)

upset(fromExpression(input), 
      nintersects = 5000, 
      nsets = 20, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 1000, 
      text.scale = c(2, 1.5, 2, 1.15, 2, 0), 
      point.size = 3.5, 
      line.size = 1,
      mainbar.y.label = "DEG intersections",
      sets.x.label = "All cell classes",
      sets.bar.color="blue"
      )
dev.off()

