import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go

# Parse arguments
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumor', type=str, help='Input tumor LR_interactions')
parser.add_argument('--normal', type=str, help='Input normal LR_interactions')

args = parser.parse_args()

print("normal file:", args.normal)
print("tumor file:", args.tumor)

# read general_info.txt in the tumor directory
info_dir = args.tumor.split("/")[0:-1]
info_file = os.path.join("/".join(info_dir), "general_info.txt")
print("Trying to read: " + info_file)

with open(info_file, "r") as f:
    for line in f:
        if "LR_resource" in line:
            resource_path = line.split(": ")[1].strip()
            break

print("Using resource file: " + resource_path)
resource = pd.read_csv(resource_path)

tumor_dir = os.path.dirname(args.tumor)
normal_dir = os.path.dirname(args.normal)

tumor_fig_dir = os.path.join(os.path.dirname(args.tumor), "figures")
normal_fig_dir = os.path.join(os.path.dirname(args.normal), "figures")

# Get name from the .h5ad file in the directory
normal_name = [f for f in os.listdir(normal_dir) if f.endswith(".h5ad")][0]
normal_name = normal_name.split(".")[0]

tumor_name = [f for f in os.listdir(tumor_dir) if f.endswith(".h5ad")][0]
tumor_name = tumor_name.split(".")[0]

print("Tumor name: " + tumor_name)
print("Normal name: " + normal_name)
print("Tumor directory: " + tumor_dir)
print("Normal directory: " + normal_dir)
print("Tumor figures will be saved to: " + tumor_fig_dir)
print("Normal figures will be saved to: " + normal_fig_dir)

print("Reading LR_interactions tables")
nlr = pd.read_csv(args.normal)
tlr = pd.read_csv(args.tumor)

# Subset to lr pairs that appear in both normal and tumor
nlr = nlr[nlr["ligand"].isin(tlr["ligand"]) & nlr["receptor"].isin(tlr["receptor"])]
tlr = tlr[tlr["ligand"].isin(nlr["ligand"]) & tlr["receptor"].isin(nlr["receptor"])]

same_module = nlr["same_module"]
nlr["consensus_module"] = "_different"
nlr.loc[same_module, "consensus_module"] = nlr.loc[same_module, "L_module"]

same_module = tlr["same_module"]
tlr["consensus_module"] = "_different"
tlr.loc[same_module, "consensus_module"] = tlr.loc[same_module, "L_module"]

print("Merging normal and tumor dataframes")
consensus_df = pd.merge(
    nlr,
    tlr,
    on=["ligand", "receptor"],
    suffixes=("_normal", "_tumor")
)[["ligand", "receptor", "consensus_module_normal", "consensus_module_tumor"]]

# Get pairs that have different modules in normal but the same in tumor
t_same_n_diff = consensus_df[
    (consensus_df["consensus_module_normal"] == "_different") & (consensus_df["consensus_module_tumor"] != "_different")
]

print('Saving interactions to file: ' + os.path.join(tumor_dir, "t_same_n_diff.csv"))
t_same_n_diff = t_same_n_diff[["ligand", "receptor"]]
t_same_n_diff.to_csv(os.path.join(tumor_dir, "t_same_n_diff.csv"), index=False)

print('Running enrichment analysis')
# Run enrichment analysis
enr = gp.enrichr(
    gene_list=t_same_n_diff["ligand"].tolist() + t_same_n_diff["receptor"].tolist(),
    gene_sets="Reactome_2022",
    background=[nlr["ligand"].tolist() + nlr["receptor"].tolist()],
    outdir=None,
)

enr= enr.results
enr = enr[enr["Adjusted P-value"] < 0.05]

enr.to_csv(os.path.join(tumor_dir, 't_same_n_diff_enrichment.csv'), index=False)

# Get pairs that have different modules in tumor but the same in normal
print("Getting pairs that have different modules in tumor but the same in normal")
n_same_t_diff = consensus_df[
    (consensus_df["consensus_module_tumor"] == "_different") & (consensus_df["consensus_module_normal"] != "_different")
]

print('Saving interactions to file: ' + os.path.join(tumor_dir, "n_same_t_diff.csv"))
n_same_t_diff = n_same_t_diff[["ligand", "receptor"]]
n_same_t_diff.to_csv(os.path.join(tumor_dir, "n_same_t_diff.csv"), index=False)

print('Running enrichment analysis')
# Run enrichment analysis
enr = gp.enrichr(
    gene_list=n_same_t_diff["ligand"].tolist() + n_same_t_diff["receptor"].tolist(),
    gene_sets="Reactome_2022",
    background=[nlr["ligand"].tolist() + nlr["receptor"].tolist()],
    outdir=None,
)

enr= enr.results
enr = enr[enr["Adjusted P-value"] < 0.05]

enr.to_csv(os.path.join(tumor_dir, 'n_same_t_diff_enrichment.csv'), index=False)

print("Getting module changes")
module_changes = pd.DataFrame(
    0,
    index=consensus_df["consensus_module_normal"].value_counts().index,
    columns=consensus_df["consensus_module_tumor"].value_counts().index,
)

for normal_module in consensus_df["consensus_module_normal"].unique():
    for tumor_module in consensus_df["consensus_module_tumor"].unique():

        module_changes.loc[normal_module, tumor_module] = sum(
            (consensus_df["consensus_module_normal"] == normal_module) & (consensus_df["consensus_module_tumor"] == tumor_module)
        )

module_changes = module_changes.fillna(0)

# Change index and columns to include type and number of pairs in each module
module_changes.index = "N" + module_changes.index.astype(str)# + ": (" + module_changes.sum(axis=1).astype(str) + " pairs)"
module_changes.columns = "T" + module_changes.columns.astype(str)# + ": (" + module_changes.sum(axis=0).astype(str) + " pairs)"

print(module_changes.head())

print("Creating Sankey plot")
# Create Sankey plot
labels = list(module_changes.index) + list(module_changes.columns)
source = []
target = []
value = []

for i, row in module_changes.iterrows():
    for j, val in row.iteritems():
        if val > 0:
            source.append(labels.index(i))
            target.append(labels.index(j))
            value.append(val)

print("Setting colors")
opacity = 0.5  # Set the fixed opacity value
t_color = (255, 0, 0) + (opacity,)  # Red
n_color = (0, 0, 255) + (opacity,)  # Blue
nan_color = (100, 100, 100) + (opacity,) # Gray

colors = []

for norm, row in module_changes.iterrows():
    for tum, val in row.iteritems():
        if val > 0:
            if (norm != "N_different") & (tum != "T_different"): # TODO: check this condition
                r, g, b, a = nan_color
                colors.append(f"rgba({r},{g},{b},{a})")
            else:
                if (norm == "N_different") & (tum == "T_different"):
                    r, g, b, a = nan_color
                    colors.append(f"rgba({r},{g},{b},{a})")
                elif (norm == "N_different"):
                    r, g, b, a = t_color
                    colors.append(f"rgba({r},{g},{b},{a})")
                elif (tum == "T_different"):
                    r, g, b, a = n_color
                    colors.append(f"rgba({r},{g},{b},{a})")
                   
r, g, b, a = t_color
tc = f"rgba({r},{g},{b},{a})"
r, g, b, a = n_color
nc = f"rgba({r},{g},{b},{a})"
node_colors = [nc] * module_changes.shape[0] + [tc] * module_changes.shape[1]
print("Generating sankey plot")

INVERT = True

if INVERT:
    source = source[::-1]
    target = target[::-1]
    value =  value[::-1]
    colors =  colors[::-1]

fig = go.Figure(
    go.Sankey(
        arrangement='snap',
        node = dict(
            pad = 15,
            thickness = 20,
            #line = dict(color = "black", width = 0.5),
            label = labels,
            color = node_colors,
            align="right",
        ),
        link = dict(
            source = source,
            target = target,
            value =  value,
            color = colors 
        )
    )
)

fig.update_layout(title_text=tumor_name, font_size=10)

print("Saving Sankey plot to " + os.path.join(tumor_fig_dir, "sankey.pdf"))
# Save figure (overwriting the dummy plot)
fig.write_image(os.path.join(tumor_fig_dir, "sankey.pdf"))

print("Done: all")
