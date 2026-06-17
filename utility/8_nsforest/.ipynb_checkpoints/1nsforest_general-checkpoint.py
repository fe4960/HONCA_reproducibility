import sys
import os
sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../nsforest/nsforesting"))
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating as ev
from nsforest import plotting as pl

data_folder=sys.argv[1]
fn=sys.argv[2]
file=data_folder + "/" + fn +  ".h5ad"
adata=sc.read_h5ad(file)

#cluster_header="subclass"
cluster_header=sys.argv[4]

#output_folder=data_folder + "/nsforest" 
output_folder=sys.argv[3] + "/nsforest_"+ cluster_header 

ns.pp.dendrogram(adata, cluster_header, save = False, output_folder = output_folder, outputfilename_suffix = cluster_header)


adata = ns.pp.prep_medians(adata, cluster_header)


adata = ns.pp.prep_binary_scores(adata, cluster_header)

plt.clf()
filename = output_folder + cluster_header + '_medians.png'
print(f"Saving median distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["medians_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{file.split("/")[-1].replace(".h5ad", "")}: {"medians_" + cluster_header} histogram')
a = plt.xlabel("medians_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')
plt.show()

plt.clf()
filename = output_folder + cluster_header + '_binary_scores.png'
print(f"Saving binary_score distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["binary_scores_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{file.split("/")[-1].replace(".h5ad", "")}: {"binary_scores_" + cluster_header} histogram')
a = plt.xlabel("binary_scores_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')
plt.show()

filename = file.replace(".h5ad", "_preprocessed.h5ad")
print(f"Saving new anndata object as...\n{filename}")
adata.write_h5ad(filename)

outputfilename_prefix = cluster_header
results = nsforesting.NSForest(adata, cluster_header, save_supplementary = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix, gene_selection="BinaryFirst_moderate")

ns.pl.boxplot(results, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

ns.pl.boxplot(results, "PPV", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

ns.pl.boxplot(results, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)


ns.pl.boxplot(results, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

ns.pl.scatter_w_clusterSize(results, "f_score", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

ns.pl.scatter_w_clusterSize(results, "PPV", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)


ns.pl.scatter_w_clusterSize(results, "recall", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)


ns.pl.scatter_w_clusterSize(results, "onTarget", save = True, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

to_plot = results.copy()

dendrogram = [] # custom dendrogram order
dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot["clusterName"] = to_plot["clusterName"].astype("category")
to_plot["clusterName"] = to_plot["clusterName"].cat.set_categories(dendrogram)
to_plot = to_plot.sort_values("clusterName")
to_plot = to_plot.rename(columns = {"NSForest_markers": "markers"})
to_plot.head()

markers_dict = dict(zip(to_plot["clusterName"], to_plot["markers"]))
markers_dict

ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = outputfilename_prefix)





