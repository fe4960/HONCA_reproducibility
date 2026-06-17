import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
# Generate random data
#data = np.random.randn(1000)  # 1000 values from a normal distribution
plt.close()
mtcars=pd.read_csv("/dfs3b/ruic20_lab/junw42/HCA_ON/data/7_DEG/genexp_sample_raw_atlas_meta_gsm_race_batch_sampleid_donor_Asian",header=0,sep="\t")
data=mtcars["age_year"]
# Create histogram
plt.figure(figsize=(10, 5))

colors = cm.get_cmap('tab20')(range(len(labels)))

#plt.hist(data, bins=50, color='lightsteelblue', edgecolor='lightsteelblue')
plt.hist(data, bins=50, color=colors[1], edgecolor=colors[1])

# Add labels and title
plt.xlabel('Age (years) ', fontsize=20)
plt.ylabel('Number of donors', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.title('Histogram Example', fontsize=16)

#plt.grid(True)
plt.tight_layout()
plt.show()

plt.savefig("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new_age_hist_python.pdf")
