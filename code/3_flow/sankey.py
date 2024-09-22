import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go

# Parse arguments
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumordir', type=str, help='Input tumordirectory')
parser.add_argument('--normaldir', type=str, help='Input normal directory')

ENRICHMENT = False

args = parser.parse_args()

tumor_dir = args.tumordir
normal_dir = args.normaldir

# Find interactions.csv files in the directories
nlr = os.path.join(normal_dir, "interactions.csv")
print("Reading normal file:", nlr)
nlr = pd.read_csv(nlr)
print(nlr.head())

tlr = os.path.join(tumor_dir, "interactions.csv")
print("Reading tumor file:", tlr)
tlr = pd.read_csv(tlr)
print(tlr.head())

# Get name from the .h5ad file in the directory
normal_name = [f for f in os.listdir(normal_dir) if f.endswith(".h5ad")][0]
normal_name = normal_name.split(".")[0]
print("Normal name: " + normal_name)

tumor_name = [f for f in os.listdir(tumor_dir) if f.endswith(".h5ad")][0]
tumor_name = tumor_name.split(".")[0]
print("Tumor name: " + tumor_name)

# Create output directories
output_dir = os.path.join("/home/lnemati/pathway_crosstalk/results/sankey", normal_name)
print("Output directory: " + output_dir)
if not os.path.exists(output_dir):
    print("Creating output directory")
    os.makedirs(output_dir)

nlr["module"] = np.where(nlr["module"].isnull(), "_different", nlr["module"])
tlr["module"] = np.where(tlr["module"].isnull(), "_different", tlr["module"])

print("Merging normal and tumor dataframes")
common_cols = ["interaction"] + [col for col in nlr.columns if col.startswith("interactor")]

consensus_df = pd.merge(
    nlr,
    tlr,
    on=common_cols,
    suffixes=("_normal", "_tumor")
)[common_cols + ["module_normal", "module_tumor"]]

print("Merged dataframe:")
print(consensus_df.head())
print("Number of interactions: " + str(consensus_df.shape[0]))

# ----- T SAME N DIFF -----

# Get pairs that have different modules in normal but the same in tumor
t_same_n_diff = consensus_df[
    (consensus_df["module_normal"] == "_different") & (consensus_df["module_tumor"] != "_different")
]

print('Saving interactions to file: ' + os.path.join(output_dir, "t_same_n_diff.csv"))
t_same_n_diff.to_csv(os.path.join(tumor_dir, "t_same_n_diff.csv"), index=False)

if ENRICHMENT:
    # Extract all genes from the columns starting with 'interactor'
    interactor_cols = [col for col in t_same_n_diff.columns if col.startswith("interactor")]

    # Get all genes from the interactors columns
    gene_list = set()
    background = set()
    for col in interactor_cols:
        gene_list = gene_list.union(set(t_same_n_diff[col].dropna()))
        background = background.union(set(consensus_df[col].dropna()))

    gene_list = list(gene_list)
    background = list(background)

    print('Running enrichment analysis')
    # Run enrichment analysis
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets="Reactome_2022",
        background=background,
        outdir=None,
    )

    enr= enr.results
    enr = enr[enr["Adjusted P-value"] < 0.05]

    enr.to_csv(os.path.join(output_dir, 't_same_n_diff_enrichment.csv'), index=False)

# ----- N SAME T DIFF -----

# Get pairs that have different modules in tumor but the same in normal
n_same_t_diff = consensus_df[
    (consensus_df["module_normal"] != "_different") & (consensus_df["module_tumor"] == "_different")
]

print('Saving interactions to file: ' + os.path.join(output_dir, "n_same_t_diff.csv"))
n_same_t_diff.to_csv(os.path.join(tumor_dir, "n_same_t_diff.csv"), index=False)

if ENRICHMENT:
    # Extract all genes from the columns starting with 'interactor'
    interactor_cols = [col for col in n_same_t_diff.columns if col.startswith("interactor")]

    # Get all genes from the interactors columns
    gene_list = set()
    background = set()
    for col in interactor_cols:
        gene_list = gene_list.union(set(n_same_t_diff[col].dropna()))
        background = background.union(set(consensus_df[col].dropna()))

    gene_list = list(gene_list)
    background = list(background)

    print('Running enrichment analysis')
    # Run enrichment analysis
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets="Reactome_2022",
        background=background,
        outdir=None,
    )

    enr= enr.results
    enr = enr[enr["Adjusted P-value"] < 0.05]

    enr.to_csv(os.path.join(output_dir, 'n_same_t_diff_enrichment.csv'), index=False)

print("Getting module changes")
module_changes = pd.DataFrame(
    0,
    index=consensus_df["module_normal"].value_counts().index,
    columns=consensus_df["module_tumor"].value_counts().index,
)

for normal_module in consensus_df["module_normal"].unique():
    for tumor_module in consensus_df["module_tumor"].unique():
        module_changes.loc[normal_module, tumor_module] = sum(
            (consensus_df["module_normal"] == normal_module) & (consensus_df["module_tumor"] == tumor_module)
        )

module_changes = module_changes.fillna(0)

# Change index and columns to include type
module_changes.index = "N" + module_changes.index.astype(str)
module_changes.columns = "T" + module_changes.columns.astype(str)

# Strip .0 if ends with it
module_changes.index = module_changes.index.str.replace("\.0", "")
module_changes.columns = module_changes.columns.str.replace("\.0", "")

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
t_color = (171, 53, 2) + (opacity,)  # sns red
n_color = (0, 114, 142) + (opacity,)  # sns blue
nan_color = (200, 202, 212) + (opacity,)  # sns grey

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

print("Saving Sankey plot to " + os.path.join(output_dir, "sankey.pdf"))
# Save figure (overwriting the dummy plot)
fig.write_image(os.path.join(output_dir, "sankey.pdf"))

print("Done: all")
