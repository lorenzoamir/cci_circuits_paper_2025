import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go
import sys

tcolor     = '#ab3502'
ncolor     = '#00728e'
graycolor  = '#4D4E4F'
graycolor2 = '#C8CAD4'

# Parse arguments
parser = argparse.ArgumentParser(description='Generate parallel categories (alluvial) plot of LR pairs')

# list of tumor and normal directories
parser.add_argument('--tumors', type=str, help='List of tumor directories separated by space', required=True)
parser.add_argument('--normals', type=str, help='List of normal directories separated by space', required=True)

args = parser.parse_args()

all_tumor_dirs = args.tumors.split(" ")
all_normal_dirs = args.normals.split(" ")

tissues_cat = []
normal_cat = []
tumor_cat = []
outcome_cat = []
counts = []
colors = []

all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')
all_interactions['same_module'] = False
all_interactions['tumor_only'] = 0 # Times the interaction is in the same module in the tumor but not in the normal
all_interactions['normal_only'] = 0 
all_interactions['both'] = 0

for tumor_dir, normal_dir in zip(all_tumor_dirs, all_normal_dirs):
    if 'testis' in normal_dir and 'testicular' in tumor_dir:
        # Testis has an extremely low modularity, so the modules are not informative
        continue

    tissue_interactions = all_interactions.copy() 
    tissue_interactions['same_module'] = False
    
    print("Tumor directory: " + tumor_dir)
    print("Normal directory: " + normal_dir)

    # Find interactions.csv files in the directories
    nlr = os.path.join(normal_dir, "interactions.csv")
    print("Reading normal file:", nlr)
    nlr = pd.read_csv(nlr)
    print(nlr.head())
    print("Number of same module interactions: " + str(sum(nlr["same_module"])))
    print()

    tlr = os.path.join(tumor_dir, "interactions.csv")
    print("Reading tumor file:", tlr)
    tlr = pd.read_csv(tlr)
    print(tlr.head())
    print("Number of same module interactions: " + str(sum(tlr["same_module"])))
    print()

    # Subset to only interactions that are in the same module in at least one of the normal or tumor networks
    print("Subsetting to only interactions that are in the same module in at least one of normal or tumor")
    keep = set(nlr[nlr['same_module']]['interaction']) | set(tlr[tlr['same_module']]['interaction'])
    tissue_interactions = tissue_interactions[tissue_interactions['interaction'].isin(keep)]
    nlr = nlr[nlr['interaction'].isin(keep)]
    tlr = tlr[tlr['interaction'].isin(keep)]
    print('Number of interactions: ' + str(tissue_interactions.shape[0]))

    # Join normal and tumor interactions with the full interaction network
    print("Merging normal and tumor interactions with the full interaction network")
    nlr = nlr.set_index('interaction')
    n_tmp = tissue_interactions.copy().set_index('interaction')
    n_tmp.loc[nlr.index, 'same_module'] = nlr['same_module']
    nlr = n_tmp
    print('Number of same module interactions: ' + str(sum(nlr['same_module'])))

    tlr = tlr.set_index('interaction') 
    t_tmp = tissue_interactions.copy().set_index('interaction')
    t_tmp.loc[tlr.index, 'same_module'] = tlr['same_module']
    tlr = t_tmp
    print('Number of same module interactions: ' + str(sum(tlr['same_module'])))
    
    # Get tissue name 
    name = normal_dir.split("/")[-3].replace("_", " ").capitalize()

    # Count occurrences of each possible outcome
    tissues_cat += [name] * 3
    normal_cat += ['Same Module', 'Same Module', 'Different Modules']
    tumor_cat += ['Same Module', 'Different Modules', 'Same Module']
    outcome_cat += ['Both', 'Normal Only', 'Tumor Only']
   
    both = tissue_interactions[nlr['same_module'].values & tlr['same_module'].values].index
    normal_only = tissue_interactions[nlr['same_module'].values & ~tlr['same_module'].values].index
    tumor_only = tissue_interactions[~nlr['same_module'].values & tlr['same_module'].values].index

    counts += [len(both), len(normal_only), len(tumor_only)]
    colors += [graycolor2, ncolor, tcolor]

    # Increment the tissue count of the interactions
    all_interactions.loc[both, 'both'] += 1
    all_interactions.loc[tumor_only, 'tumor_only'] += 1
    all_interactions.loc[normal_only, 'normal_only'] += 1
    

# Save categories and counts as a csv file
df = pd.DataFrame({'Tissue': tissues_cat, 'Normal': normal_cat, 'Tumor': tumor_cat, 'Outcome': outcome_cat, 'Counts': counts})
df.to_csv('/home/lnemati/pathway_crosstalk/results/flow/flow_diagram_data.csv', index=False)


# Define index that measures the rewiring of the interaction
# T - N / T + N + B
sums = all_interactions['tumor_only'] + all_interactions['normal_only'] + all_interactions['both']
all_interactions = all_interactions[sums > 0]
sums = all_interactions['tumor_only'] + all_interactions['normal_only'] + all_interactions['both']
diff = all_interactions['tumor_only'] - all_interactions['normal_only']

# Remove interactions that never appear
all_interactions['rewiring_index'] = diff / sums

# Save the interactions
all_interactions.sort_values(by=['rewiring_index'], ascending=False, inplace=True)
all_interactions.to_csv('/home/lnemati/pathway_crosstalk/results/flow/interactions_with_counts.csv', index=False)

fs = 40 # font size

# Make parallel categories plot, use colors to distinguish between normal and tumor
fig = go.Figure(data=[go.Parcats(
    dimensions=[
        {'label': 'Tissue',
         'values': tissues_cat},
        {'label': 'Normal',
         'values': normal_cat},
        {'label': 'Tumor',
         'values': tumor_cat},
        {'label': 'Outcome',
         'values': outcome_cat}
    ],
    counts=counts,
    line={'shape': 'hspline', 'color': colors},
    )]
)

fig.update_layout(
    #title_text='Flow of LR pairs',
    margin=dict(l=200, r=200, t=50, b=50),
    font_size=fs)

# Save figure
output_dir = '/home/lnemati/pathway_crosstalk/results/flow'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save figure make the figure big
fig.write_image(os.path.join(output_dir, "flow.pdf"), width=2400, height=1400)
fig.write_image(os.path.join(output_dir, "flow.png"), width=2400, height=1400)

print("Done: flow.py")
