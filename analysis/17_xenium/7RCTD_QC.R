resultsdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_subclass1_01"
results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                     results$results_df)
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)


resultsdir="/dfs3b/ruic20_lab/junw42/HCA_ON/data/8_spatial/Xenium/20250620_PP_comb/ONONH_comb/spacexr/oligo_subclass_015"
results <- RCTD@results
# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- RCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- RCTD@spatialRNA
#resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)
plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                     results$results_df)
plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
