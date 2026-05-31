#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import pandas as pd
from scenicplus.RSS import plot_rss
import matplotlib.backends.backend_pdf

indata=pd.read_csv(infile, sep='\t', header=0, index_col=indexcol)
print(f"{indata=}", flush=True)

plot_rss(
	data_matrix=indata,
	top_n=ntop,
	num_columns=numcolumn,
	figsize=(width, height),
	save=f"{outdir}/{bname}.pdf",
)
