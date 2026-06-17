import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm
plt.close("all")

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
cellnum=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new.csv",header=0,index_col=0)

cellnum1=cellnum.loc[-cellnum["major class"].isin(["BC","Rod","Cone","HC","AC","RGC","Adipocyte","Pigmented_cell","RPE","Schwann_cell"]),]

# Sample data
#df <- data.frame(
#  category = c("A", "B", "C", "D"),
#  value = c(10, 15, 100, 12)  # 'C' is an outlier
#)
#cellnum1=cellnum1[order(cellnum1$cell.number,decreasing = TRUE),]
#######cellnum1=cellnum1.sort_values(by="cell.number",ascending=False)
# Sample data
categories = cellnum1["major class"]  #['A', 'B', 'C', 'D']
values = cellnum1["cell number"] #[10, 15, 100, 12]  # 'C' is a large outlier

# Create the figure and two subplots (stacked vertically)
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 6), gridspec_kw={'height_ratios': [1, 2]})

# Plot on both axes
#ax1.bar(categories, values, color='skyblue')
#ax2.bar(categories, values, color='skyblue')
colors = cm.get_cmap('tab20')(range(len(categories)))
ax1.bar(categories, values, color=colors)
ax2.bar(categories, values, color=colors)

# Set y-axis limits to create the break
ax1.set_ylim(100000, 280000)   # Upper part, showing only the top of the tall bar
ax2.set_ylim(0, 50000)     # Lower part, showing the shorter bars

# Hide the spines between the two plots
ax1.spines['bottom'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax1.tick_params(bottom=False)

# Add diagonal lines to indicate axis break
d = 0.015  # Size of the diagonal lines
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((-d, +d), (-d, +d), **kwargs)        # Top-left diagonal
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # Top-right diagonal

kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)        # Bottom-left diagonal
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # Bottom-right diagonal

# Add labels and title
plt.xlabel("Major class", fontsize=15) #,rotation=90, labelpad=10)
ax2.set_ylabel("Cell number",  fontsize=15)
fig.suptitle("Truncated Bar Plot (Broken Axis)")
plt.setp(ax2.get_xticklabels(), rotation=90, ha='center', fontsize=15)

plt.tight_layout()
plt.show()
fig.savefig("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new_truncated_bar_plot.pdf", format='pdf', bbox_inches='tight')
