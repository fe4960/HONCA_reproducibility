#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import mudata
from scenicplus.plotting.dotplot import heatmap_dotplot
import pandas as pd

if infile.endswith('.h5mu'):
	indata=mudata.read(infile)
else:
	print(f"Error: format is not supported for {infile}")
	sys.exit(-1)

print(f"{indata=}", flush=True)

# Subset features
featurenames=pd.read_csv(namefile, sep='\t', header=None).iloc[:, 0].tolist() if namefile else None
print(f"{featurenames=}", flush=True)

plot=heatmap_dotplot(
	scplus_mudata=indata,
	color_modality="direct_gene_based_AUC",
	size_modality="direct_region_based_AUC",
	group_variable=f"scRNA_counts:{label}",
	eRegulon_metadata_key="direct_e_regulon_metadata",
	color_feature_key="Gene_signature_name",
	size_feature_key="Region_signature_name",
	feature_name_key="eRegulon_name",
	sort_data_by="direct_gene_based_AUC",
	orientation="horizontal",
	figsize=(width, height),
	split_repressor_activator=True,
	save=None,
	subset_feature_names=featurenames,
	# save=f"{outdir}/{bname}.pdf", # limited by figure size
)
plot.save(f"{outdir}/{bname}.{format}", width=width, height=height, limitsize=False, format=format) # bug fix figsize
