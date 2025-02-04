import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go

tcolor     = '#ab3502'
ncolor     = '#00728e'
graycolor2 = '#C8CAD4'

def hex_to_rgba(hex_color, alpha=0.7):
    """
    Convert a hex color to an RGBA string with specified transparency.

    Parameters:
        hex_color (str): The hex color string (e.g., '#ab3502').
        alpha (float): Transparency level (0.0 to 1.0, default is 0.5).

    Returns:
        str: The RGBA color string (e.g., 'rgba(171, 53, 2, 0.5)').
    """
    # Remove '#' if present
    hex_color = hex_color.lstrip('#')
    
    # Convert hex to RGB components
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
    
    # Return the RGBA string
    return f'rgba({r}, {g}, {b}, {alpha})'

tcolor= hex_to_rgba(tcolor, 0.5)
ncolor= hex_to_rgba(ncolor, 0.5)
graycolor2= hex_to_rgba(graycolor2, 0.5)

# Parse arguments
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumor', type=str, help='Input tumor LR_interactions')
parser.add_argument('--normal', type=str, help='Input normal LR_interactions')

args = parser.parse_args()

print("normal file:", args.normal)
print("tumor file:", args.tumor)

tumor_dir = os.path.dirname(os.path.dirname(args.tumor))
normal_dir = os.path.dirname(os.path.dirname(args.normal))
print("tumor_dir:", tumor_dir)
print("normal_dir:", normal_dir)

tumor_fig_dir = os.path.join(tumor_dir, "figures")
normal_fig_dir = os.path.join(normal_dir, "figures")
print("tumor_fig_dir:", tumor_fig_dir)
print("normal_fig_dir:", normal_fig_dir)

tumor_module_info = os.path.join(tumor_dir, "module_info.csv")
normal_module_info = os.path.join(normal_dir, "module_info.csv")
print("tumor_module_info:", tumor_module_info)
print("normal_module_info:", normal_module_info)
tumor_module_info = pd.read_csv(tumor_module_info)
normal_module_info = pd.read_csv(normal_module_info)

# Get name from the .h5ad file in the directory
normal_name = [f for f in os.listdir(normal_dir) if f.endswith(".h5ad")][0]
normal_name = normal_name.split(".")[0]
print("normal_name:", normal_name)

tumor_name = [f for f in os.listdir(tumor_dir) if f.endswith(".h5ad")][0]
tumor_name = tumor_name.split(".")[0]
print("tumor_name:", tumor_name)

print("Tumor name: " + tumor_name)
print("Normal name: " + normal_name)
print("Tumor directory: " + tumor_dir)
print("Normal directory: " + normal_dir)
print("Tumor figures will be saved to: " + tumor_fig_dir)
print("Normal figures will be saved to: " + normal_fig_dir)

print("Reading interactions")
nlr = pd.read_csv(args.normal)
tlr = pd.read_csv(args.tumor)

# Set nans to -1
n_nan = nlr['module'].isnull()
t_nan = tlr['module'].isnull()
nlr['module'] = nlr['module'].fillna(-1)
tlr['module'] = tlr['module'].fillna(-1)

# Make sure the module columns are integers
nlr["module"] = nlr["module"].astype(int)
tlr["module"] = tlr["module"].astype(int)
print("Normal modules value counts:")
print(nlr["module"].value_counts())
print("Tumor modules value counts:")
print(tlr["module"].value_counts())

# Use module_info to convert module names to their rank
print("Converting module names to their rank")
tmodule2rank = dict(zip(tumor_module_info["module"], tumor_module_info["rank"]))
nmodule2rank = dict(zip(normal_module_info["module"], normal_module_info["rank"]))
nlr["module"] = nlr["module"].map(nmodule2rank)
tlr["module"] = tlr["module"].map(tmodule2rank)

# Write 'Different<br>Modules' if the modules are different (nans)
nlr["module"] = np.where(n_nan, "Different<br>Normal<br>Modules", nlr["module"]).astype(str)
tlr["module"] = np.where(t_nan, "Different<br>Tumor<br>Modules", tlr["module"]).astype(str)
# Strip .0 if ends with it
nlr["module"] = nlr["module"].str.replace("\.0", "")
tlr["module"] = tlr["module"].str.replace("\.0", "")
# Print value counts
print("Normal modules value counts:")
print(nlr["module"].value_counts())
print("Tumor modules value counts:")
print(tlr["module"].value_counts())

print("Merging normal and tumor dataframes")
consensus_df = pd.merge(
    nlr,
    tlr,
    on='interaction',
    suffixes=("_normal", "_tumor")
)[['interaction', 'module_normal', 'module_tumor']]

print("Merged dataframe:")
print(consensus_df.head())
print("Number of interactions: " + str(consensus_df.shape[0]))

# Get pairs that have different modules in normal but the same in tumor
t_same_n_diff = consensus_df[
    (consensus_df["module_normal"] == "Different<br>Normal<br>Modules") & (consensus_df["module_tumor"] != "Different<br>Tumor<br>Modules")
]

print('Saving interactions to file: ' + os.path.join(tumor_dir, "t_same_n_diff.csv"))
t_same_n_diff.to_csv(os.path.join(tumor_dir, "t_same_n_diff.csv"), index=False)

# Get pairs that have different modules in tumor but the same in normal
n_same_t_diff = consensus_df[
    (consensus_df["module_normal"] != "Different<br>Normal<br>Modules") & (consensus_df["module_tumor"] == "Different<br>Tumor<br>Modules")
]
print('Saving interactions to file: ' + os.path.join(tumor_dir, "n_same_t_diff.csv"))
n_same_t_diff.to_csv(os.path.join(tumor_dir, "n_same_t_diff.csv"), index=False)

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

print('Module changes:')
print(module_changes.head())

# Strip .0 if ends with it
print("Stripping .0 from module names")
module_changes.index = module_changes.index.str.replace(".0", "")
module_changes.columns = module_changes.columns.str.replace(".0", "")

# Reorder index and columns so that "Different<br>Modules" is the first columns, the rest are sorted by their rank
# NOT ALPHABETICALLY, so 1,2,3,4... not 10,11,12...
print('Sorting modules')
# Convert to float for sorting, but then go back to strings
t_sorted = sorted([int(col) for col in module_changes.columns if col != "Different<br>Tumor<br>Modules"])
t_sorted = [str(int(col)) for col in t_sorted]
module_changes = module_changes[
    ["Different<br>Tumor<br>Modules"] + t_sorted
]

n_sorted = sorted([float(idx) for idx in module_changes.index if idx != "Different<br>Normal<br>Modules"])
n_sorted = [str(int(idx)) for idx in n_sorted]
module_changes = module_changes.loc[
    ["Different<br>Normal<br>Modules"] + n_sorted
]
# Add N and T to the index and columns (only if they are not "Different<br>Modules")
module_changes.columns = np.where(
    module_changes.columns == "Different<br>Tumor<br>Modules",
    "Different<br>Tumor<br>Modules",
    "T" + module_changes.columns.astype(str)
)

module_changes.index = np.where(
    module_changes.index == "Different<br>Normal<br>Modules",
    "Different<br>Normal<br>Modules",
    "N" + module_changes.index.astype(str)
)

print(module_changes.head())

print("Creating Sankey plot")
# Create Sankey plot
labels = list(module_changes.index) + list(module_changes.columns)
source = []
target = []
value = []

# Set to 0 the cell corresponding to both in different modules
print('Discardin different modules in both conditions')
module_changes.loc["Different<br>Normal<br>Modules", "Different<br>Tumor<br>Modules"] = 0

for i, row in module_changes.iterrows():
    for j, val in row.items():
        source.append(labels.index(i))
        target.append(labels.index(j))
        value.append(val)

print("Setting colors")
colors = []

for norm, row in module_changes.iterrows():
    for tum, val in row.items():
        if (norm != "Different<br>Normal<br>Modules") & (tum != "Different<br>Tumor<br>Modules"):
            #r, g, b, a = nan_color
            #colors.append(f"rgba({r},{g},{b},{a})")
            colors.append(graycolor2)
        else:
            if (norm == "Different<br>Normal<br>Modules") & (tum == "Different<br>Tumor<br>Modules"):
                #r, g, b, a = nan_color
                #colors.append(f"rgba({r},{g},{b},{a})")
                colors.append(graycolor2)	
            elif (norm == "Different<br>Normal<br>Modules"):
                #r, g, b, a = t_color
                #colors.append(f"rgba({r},{g},{b},{a})")
                colors.append(tcolor)
            elif (tum == "Different<br>Tumor<br>Modules"):
                #r, g, b, a = n_color
                #colors.append(f"rgba({r},{g},{b},{a})")
                colors.append(ncolor)
                   
# Different modules are gray, the rest are normal or tumor colors
node_colors = [graycolor2] + \
              [ncolor] * (module_changes.shape[0] - 1) + \
              [graycolor2] + \
              [tcolor] * (module_changes.shape[1] - 1)

# Defining positions

# Normal nodes on the left, tumor nodes on the right
x = [0.01] * module_changes.shape[0] + [0.99] * module_changes.shape[1]

# Order nodes vertically by rank
y = list(np.linspace(0.01, 0.99, module_changes.shape[0])) + list(np.linspace(0.01, 0.99, module_changes.shape[1])) 

print("Generating sankey plot")

# print labels and positions
print('x:', x)
print([(posx, lab) for posx, lab in zip(x, labels)])
print('y:', y)
print([(posy, lab) for posy, lab in zip(y, labels)])

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
            pad = 20,
            thickness = 20,
            line = dict(color = "black", width = 1.5),
            label = labels,
            color = node_colors,
            x = x,
            y = y,
            #align="right",
        ),
        link = dict(
            source = source,
            target = target,
            value =  value,
            color = colors, 
            line=dict(color="black", width=0),
        ),
    )
)

title=tumor_name.replace("_", " ").title()
fig.update_layout(
    title=dict(text=title),
    font=dict(size=16, color='black')
)

print("Saving Sankey plot to " + os.path.join(tumor_fig_dir, "sankey.pdf"))
# Save figure (overwriting the dummy plot)
fig.write_image(os.path.join(tumor_fig_dir, "sankey.pdf"))
fig.write_image(os.path.join(tumor_fig_dir, "sankey.svg"))
fig.write_image(os.path.join(tumor_fig_dir, "sankey.png"))

print("Done: sankey.py")
