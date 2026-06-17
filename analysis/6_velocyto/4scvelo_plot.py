import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib.pyplot as plt
import sys
#scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

indir=sys.argv[1]
h5ad=sys.argv[2]
dir1=f"{indir}/veloc"

adata=sc.read(f'{dir1}/{h5ad}_velocity.h5ad')
scv.pl.velocity_embedding_stream(adata, color="subclass",basis='umap', save=f'{dir1}/{h5ad}_velocity.png')

