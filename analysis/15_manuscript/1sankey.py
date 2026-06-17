import plotly.graph_objects as go
import pandas as pd
data=pd.read_csv("HCA_ON/scripts/15_manuscript/list/sankey_data",header=None,sep="\t")
label=pd.read_csv("HCA_ON/scripts/15_manuscript/list/sankey_label",header=None,sep="\t")

fig = go.Figure(data=[go.Sankey(
    node=dict(
        pad=40,
        thickness=40,
        line=dict(color="black", width=0.5),
        label=label[0],
#        palette="tab20"
#        color="lightblue"
    ),
    link=dict(
        source=data[0],
        target=data[1],
        value=data[2]
    )
)])

#fig.update_layout(title_text="Basic Sankey Diagram", font_size=20,width=1000, height=600 )
fig.update_layout( font=dict(size=40, color="black",family="Arial", weight=500),width=1700, height=2100 )
#fig.update_layout( font_family="Arial", font_size=30, font_color="black",width=1200, height=1500 )

# Save to PNG file
fig.write_image("HCA_ON/data/15_manuscript/sankey_diagram.pdf") #, width=1000, height=1000)
