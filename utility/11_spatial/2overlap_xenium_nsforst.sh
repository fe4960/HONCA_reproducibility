cut -f 1 HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata | sort -u > HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene
awk -F "," '{print $2}' /dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest/nsforest_celltypecelltype_markers.csv | sort -u > HCA_ON/data/8_spatial/Xenium/nsforest_celltypecelltype_markers_gene

dif.pl HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene HCA_ON/data/8_spatial/Xenium/nsforest_celltypecelltype_markers_gene HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene_u_ns_mk  HCA_ON/data/8_spatial/Xenium/nsforest_celltypecelltype_markers_gene_u_xenium_prob HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene_overlap_ns_mk

grep -f  HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene_overlap_ns_mk  /dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/nsforest/nsforest_celltypecelltype_markers.csv > HCA_ON/data/8_spatial/Xenium/Xenium_hSkin_v1_metadata_gene_overlap_ns_mk_detail
