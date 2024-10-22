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
parser.add_argument('--tissues', type=str, help='List of tissue directories separated by space', required=True)

args = parser.parse_args()

# Inside each tissue directory, there are two directories: normal and tumor. Find all subdirectories of normal and tumor
print('Input tissues:', args.tissues)

tumor_dirs = {}
normal_dirs = {}

all_tissues = args.tissues.split(" ")

# Remove empty strings and trailing characters
all_tissues = [x.strip() for x in all_tissues if x.strip() != '']

for tissue_dir in all_tissues:
    # Get tissue name from path
    tissue_name = tissue_dir.split("/")[-1]

    # Get tumor subdirectories
    tumor_dirs[tissue_name] = [os.path.join(tissue_dir, 'tumor', x) for x in os.listdir(os.path.join(tissue_dir, 'tumor'))]
    # Get normal subdirectories
    normal_dirs[tissue_name] = [os.path.join(tissue_dir, 'normal', x) for x in os.listdir(os.path.join(tissue_dir, 'normal'))]

tissues_cat = []
normal_cat = []
tumor_cat = []
outcome_cat = []
counts = []
colors = []

all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv', index_col='interaction')
all_interactions['both_score'] = 0
all_interactions['normal_score'] = 0
all_interactions['tumor_score'] = 0

def format_string(string, newline=False):
    '''
    replaces _ with space and capitalizes the first letter of each word.
    finds which of the spaces can be replaced with a newline as to make the
    two resulting lines of text as close to equal length as possible.
    '''
    # Replace _ with space
    string = string.replace('_', ' ').title()
    
    if newline:
        # If string has less than 20 characters, don't bother
        if len(string) < 10:
            return string
        # Find the space that makes the two lines of text as close to equal length as possible
        diffs = []

        for i, words in enumerate(string.split()):
            diffs.append(abs(len(' '.join(string.split()[:i])) - len(' '.join(string.split()[i:]))))
            
        # Find which space corresponds to the minimum difference and replace it with <br>
        n_space = diffs.index(min(diffs)) # e.g. the second ' ' in the string
        string = ' '.join(string.split()[:n_space]) + '<br>' + ' '.join(string.split()[n_space:])

    return string

#for tumor_dir, normal_dir in zip(all_tumor_dirs, all_normal_dirs):
for tissue_dir in all_tissues:
    # Testis has an extremely low modularity, so the modules are not informative
    if 'testis' in tissue_dir:
        continue

    tissue_name = tissue_dir.split("/")[-1] 
    print(tissue_name)
    
    tissue_interactions = all_interactions.copy() 
    tissue_interactions['same_module'] = False

    # Find all interactions.csv files in the normal and tumor directories
    normal_dfs = [os.path.join(n_dir, 'interactions.csv') for n_dir in normal_dirs[tissue_name]]
    tumor_dfs = [os.path.join(t_dir, 'interactions.csv') for t_dir in tumor_dirs[tissue_name]]

    print('Normal tissues:', len(normal_dfs))
    print('Tumor tissues:', len(tumor_dfs))
    
    # Read the dataframes
    normal_dfs = [pd.read_csv(df, index_col='interaction') for df in normal_dfs]
    tumor_dfs = [pd.read_csv(df, index_col='interaction') for df in tumor_dfs]

    # Print shapes
    print('Normal shapes:', [df.shape for df in normal_dfs])
    print('Tumor shapes:', [df.shape for df in tumor_dfs])
   
    print('Normal same module:', [sum(df['same_module']) for df in normal_dfs])
    print('Tumor same module:', [sum(df['same_module']) for df in tumor_dfs])

    ## Only keep interactions that are in same module in at least one of the normal or tumor networks
    #keep = set()

    #for df in normal_dfs + tumor_dfs:
    #    keep |= set(df[df['same_module']].index)
   
    #print('Found', len(keep), 'interactions to keep')
    #print('Subsetting')
    ## Subset the tissue_interactions
    #tissue_interactions = tissue_interactions.loc[keep]

    # Merge normal and tumor interactions with the full interaction network
    print('Merging')
    new_dfs = []
    for df in normal_dfs:
        n_tmp = tissue_interactions.copy()
        n_tmp.loc[df.index, 'same_module'] = df['same_module']
        new_dfs.append(n_tmp)
    normal_dfs = new_dfs

    new_dfs = []
    for df in tumor_dfs:
        t_tmp = tissue_interactions.copy()
        t_tmp.loc[df.index, 'same_module'] = df['same_module']
        new_dfs.append(t_tmp)
    tumor_dfs = new_dfs
   
    # Make sure all dfs follow the same order
    normal_dfs = [df.loc[tissue_interactions.index] for df in normal_dfs]
    tumor_dfs = [df.loc[tissue_interactions.index] for df in tumor_dfs]

    print('Normal shapes:', [df.shape for df in normal_dfs])
    print('Tumor shapes:', [df.shape for df in tumor_dfs])

    print('Normal same module:', [sum(df['same_module']) for df in normal_dfs])
    print('Tumor same module:', [sum(df['same_module']) for df in tumor_dfs])

    # Sum the number of times the interaction is in the same module in the normal and tumor networks
    tissue_interactions['tot_normal'] = sum([df['same_module'] for df in normal_dfs])
    tissue_interactions['tot_tumor'] = sum([df['same_module'] for df in tumor_dfs])

    # Also get the fraction of times the interaction is in the same module
    tissue_interactions['frac_normal'] = tissue_interactions['tot_normal'] / len(normal_dfs)
    tissue_interactions['frac_tumor'] = tissue_interactions['tot_tumor'] / len(tumor_dfs)
    
    print('First 5 interactions:')
    print('Aggregated:')
    print(tissue_interactions[['tot_normal', 'tot_tumor', 'frac_normal', 'frac_tumor']].head())
    
    print('Normal:')
    print([normal_df['same_module'].head() for normal_df in normal_dfs])

    print('Tumor:')
    print([tumor_df['same_module'].head() for tumor_df in tumor_dfs])

    # Get flow for each outcome
    tissues_cat += [format_string(tissue_name, newline=False)] * 3
    normal_cat += ['Same<br>Module', 'Same<br>Module', 'Different<br>Modules']
    tumor_cat += ['Same<br>Module', 'Different<br>Modules', 'Same<br>Module']
    outcome_cat += ['Both', 'Normal<br>Only', 'Tumor<br>Only']
    
    # Both is the fraction of tissues in which the interaction is in the same module in both normal and tumor
    # So its the minimum of the two fractions
    both = tissue_interactions[['frac_normal', 'frac_tumor']].min(axis=1)
    normal_only = tissue_interactions['frac_normal'] - both
    tumor_only = tissue_interactions['frac_tumor'] - both
    
    # Save fractions to all_interactions
    all_interactions.loc[tissue_interactions.index, 'both_score'] += both
    all_interactions.loc[tissue_interactions.index, 'normal_score'] += normal_only
    all_interactions.loc[tissue_interactions.index, 'tumor_score'] += tumor_only

    # Sum values
    both = both.sum()
    normal_only = normal_only.sum()
    tumor_only = tumor_only.sum()

    counts += [both, normal_only, tumor_only]
    colors += [graycolor2, ncolor, tcolor]

n_tissues = len(tissues_cat) // 3

# Save categories and counts as a csv file
df = pd.DataFrame({'Tissue': tissues_cat, 'Normal': normal_cat, 'Tumor': tumor_cat, 'Outcome': outcome_cat, 'Counts': counts})
df.to_csv('/home/lnemati/pathway_crosstalk/results/flow/flow_diagram_data.csv', index=False)

# Define index that measures the rewiring of the interaction
#sums = all_interactions['tumor_only'] + all_interactions['normal_only'] + all_interactions['both']
#all_interactions = all_interactions[sums > 0]

# Rewiring index is T**2 - N**2 / N_tissues**2
# 1 if all tumor, -1 if all normal, 0 if equal
# also close to 0 if T and N are smaller than N_tissues
#all_interactions['rewiring_index'] = (all_interactions['tumor_only']**2 - all_interactions['normal_only']**2) / n_tissues**2

# Save the interactions
#all_interactions = all_interactions.sort_values('rewiring_index', ascending=False)

# Remove interactions that only have 0 values
keep = all_interactions[(all_interactions['both_score'] > 0) | (all_interactions['normal_score'] > 0) | (all_interactions['tumor_score'] > 0)].index
all_interactions = all_interactions.loc[keep]
all_interactions['diff'] = all_interactions['tumor_score'] - all_interactions['normal_score']

# Sort by difference between tumor and normal
all_interactions = all_interactions.sort_values('diff', ascending=False)
# Only save scores
all_interactions = all_interactions[['both_score', 'normal_score', 'tumor_score']]
all_interactions.to_csv('/home/lnemati/pathway_crosstalk/results/flow/interactions_with_counts.csv')

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

# Larger margins, larger font, black font color
fig.update_layout(
    #title_text='Flow of LR pairs',
    margin=dict(l=250, r=120, t=50, b=50),
    font_size=45,
    font_color='black'
)

# Save figure
output_dir = '/home/lnemati/pathway_crosstalk/results/flow'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Save figure make the figure big
fig.write_image(os.path.join(output_dir, "flow.pdf"), width=3000, height=1800)
fig.write_image(os.path.join(output_dir, "flow.png"), width=3000, height=1800)

print("Done: flow.py")
