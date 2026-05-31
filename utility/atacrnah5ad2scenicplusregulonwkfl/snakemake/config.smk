import yaml
import json # for debug

# Parameters
(
	atac,
	rna,
	bname,
	nowtimestr,
	)=(
		config['atac'],
		config['rna'],
		config['bname'],
		config['nowtimestr'],
		)

# Default parameters
default_params={
	'scrnah5adsubsetsamplingbykey': {
		'condaenv': None,
		'atackey': 'celltype',
		'rnakey': 'celltype',
		'ncell': 2000,
		'seed': 42,
	},
	'scrnah5adduplicateobs': {
		'condaenv': None,
		'key': '_scenicplus_celltype',
	},
	'scrnah5adobscolumn2strreplace': {
		'condaenv': None,
		'pattern': '_',
		'replace': '-',
	},
	'atacpeakset2cistopicmodel': {
		'condaenv': None,
		'blacklist': 'hg38',
		'numthreads': 12,
		'mallet': '/dfs3b/ruic20_lab/jinl14/.sync/tools/Mallet/Mallet-202108/bin/mallet',
		'memory': 200,
		'ntopics': [2] + list(range(5, 151, 5)),
		'niter': 500,
		'seed': 42,
		'alpha': 50,
		'eta': 0.1,
		'filterpeakset': True,
	},
	'ataccombinecistopicmodel2evaluate': {
		'condaenv': None,
	},
	'ataccistopicobject2addldamodel': {
		'condaenv': None,
	},
	'ataccistopicobject2celltopicheatmap': {
		'condaenv': None,
		'width': 10,
		'height': 8,
		'format': 'png',
	},
	'ataccistopicobject2runumap': {
		'condaenv': None,
		'width': 8,
		'height': 8,
		'format': 'png',
		'label': ['Sample'],
	},
	'ataccistopicobject2regionsets': {
		'condaenv': None,
		'ntop': 3000,
		'width': 20,
		'height': 20,
		'format': 'png',
		'qvalue': 0.05,
		'log2foldchange': 1.5,
		'numthreads': 12,
		'rundar': True,
	},
	'ataccistopic2runscenicpluspipeline': {
		'condaenv': None,
		'numthreads': 24,
		'ctxdb': '/dfs3b/ruic20_lab/jinl14/resource/cisTarget/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather',
		'demdb': '/dfs3b/ruic20_lab/jinl14/resource/cisTarget/hg38_screen_v10_clust.regions_vs_motifs.scores.feather',
		'motif': '/dfs3b/ruic20_lab/jinl14/resource/cisTarget/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl',
		'genome': 'hg38',
		'dem_qthr': 0.05,
		'dem_log2fcthr': 1.0,
	},
	'mudatah5mu2unstxt': {
		'condaenv': None,
	},
	'scenicplush5mu2eregulonumap': {
		'condaenv': None,
		'width': 6,
		'height': 6,
		'seed': 42,
		'format': 'png',
		'notitle': False,
	},
	'scenicplush5mu2eregulonrss': {
		'condaenv': None,
	},
	'scenicplush5mu2regulonqcpseudobulks': {
		'condaenv': None,
		'modkey_gene': 'direct_gene_based_AUC',
		'modkey_region': 'direct_region_based_AUC',
		'nsample': 10,
		'nbulk': 100,
		'seed': 42,
	},
	'scenicplusunstxt2generegionfeaturemap': {
		'condaenv': None,
		'src': 'Gene_signature_name',
		'target': 'Region_signature_name',
	},
	'scenicplusregulonpseudobulks2correlation': {
		'condaenv': None,
	},
	'scenicplusreguloncorrelation2subsetfeature': {
		'condaenv': None,
		'regulonname': 'eRegulon_name',
		'qthr': 0.05,
		'rhothr': 0.5,
		'region_count': 20,
	},
	'scenicplush5mu2eregulonheatmap': {
		'condaenv': None,
		'width': 36,
		'height': 12,
	},
	'scenicplusrsstxt2plotrss': {
		'condaenv': None,
		'ntop': 10,
		'width': 3.5,
		'height': 4,
		'numcolumn': 4,
	},
}

# Assign default params to config if not existing
for key, paramdict in default_params.items():
	if key not in config:
		config[key]=paramdict
	else:
		for subkey, value in paramdict.items():
			if subkey not in config[key]:
				config[key][subkey]=value
			else:
				print(f"Info: {key} and {subkey} are in config. So, skipping default assignment.", flush=True)

# Overwrite condaenv
if config['condaenv']=='None':
	config['condaenv']=None
if config['condaenv']:
	for key, paramdict in config.items():
		if isinstance(paramdict, dict) and 'condaenv' in paramdict:
			paramdict['condaenv']=config['condaenv']

print(json.dumps(config, indent=4), flush=True)

# save parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
