import pandas as pd
import numpy as np
import sys
outdir=sys.argv[2]
bname=sys.argv[1]
conf_mat=pd.read_csv(f'{outdir}/{bname}_conf_mat_table.csv',header=0, index_col=0)

#conf_mat=pd.read_csv(f'{outdir}/{bname}_sample_conf_mat_table.csv',header=0, index_col=0)
row_max=conf_mat.max(axis=1)
dominate_name=row_max[row_max > 0.8].index
print(dominate_name)
conf_mat.idxmax(axis=1).to_csv(f"{outdir}/{bname}_conf_mat_table_max.csv")

#dominate_name.to_csv(f'{outdir}/{bname}_dominate_sample.csv',header=True, index=True)
