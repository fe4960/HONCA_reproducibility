import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams["font.sans-serif"] = ["Nimbus Sans"]
# Sample data
colors = cm.get_cmap('tab20')(range(len(labels)))

plt.close()

labels = ['PMID:37566633','Unpublished']
sizes = [142708, 739477]

colors = cm.get_cmap('tab10')(range(len(labels)))

#colors = ['gold', 'skyblue', 'lightcoral', 'limegreen']

# Custom function to show both label and value
def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return f'{pct:.1f}%\n({val})'
    return my_autopct

# Create pie chart
plt.figure(figsize=(6, 6))
plt.pie(
    sizes, 
    labels=labels, 
    colors=colors, 
    autopct=make_autopct(sizes), 
    startangle=90,
    textprops={'fontsize': 20}
)
plt.title("Study distribution", fontsize=20)
plt.axis('equal')  # Keeps the pie chart circular

plt.show()
plt.savefig("/dfs3b/ruic20_lab/junw42/HCA_ON/data/5_refine_major/scvi/major/figures/major_ON_ONH_clean_new_study_piechart_python.pdf")
