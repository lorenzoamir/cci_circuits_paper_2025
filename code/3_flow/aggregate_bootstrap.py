import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='Merge together the results of the bootstrap same_module analysis')

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

all_interactions = pd.read_csv('/home/lnemati/pathway_crosstalk/results/flow/bootstrap/bootstrap_all.csv', index_col=0)
all_interactions['both_score'] = 0
all_interactions['normal_score'] = 0
all_interactions['tumor_score'] = 0


#for tumor_dir, normal_dir in zip(all_tumor_dirs, all_normal_dirs):
for tissue_dir in all_tissues:
    # Testis has an extremely low modularity, so the modules are not informative
    if 'testis' in tissue_dir:
        continue

    tissue_name = tissue_dir.split("/")[-1] 
    print(tissue_name)
    
    tissue_interactions = all_interactions.copy() 
    tissue_interactions['same_module'] = False

    # Find all interactions_bootstrap.csv files in the normal and tumor directories
    normal_dfs = [os.path.join(n_dir, 'interactions_bootstrap.csv') for n_dir in normal_dirs[tissue_name]]
    tumor_dfs = [os.path.join(t_dir, 'interactions_bootstrap.csv') for t_dir in tumor_dirs[tissue_name]]

    print('Normal tissues:', len(normal_dfs))
    print('Tumor tissues:', len(tumor_dfs))
    
    # Read the dataframes
    normal_dfs = [pd.read_csv(df, index_col=0) for df in normal_dfs]
    tumor_dfs = [pd.read_csv(df, index_col=0) for df in tumor_dfs]

    # Print shapes
    print('Normal shapes:', [df.shape for df in normal_dfs])
    print('Tumor shapes:', [df.shape for df in tumor_dfs])
   
    print('Normal same module:', [sum(df['same_module']) for df in normal_dfs])
    print('Tumor same module:', [sum(df['same_module']) for df in tumor_dfs])

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

    # Both is the fraction of tissues in which the interaction is in the same module in both normal and tumor
    # So its the minimum of the two fractions
    both = tissue_interactions[['frac_normal', 'frac_tumor']].min(axis=1)
    normal_only = tissue_interactions['frac_normal'] - both
    tumor_only = tissue_interactions['frac_tumor'] - both
    
    # Save fractions to all_interactions
    all_interactions.loc[tissue_interactions.index, 'both_score'] += both
    all_interactions.loc[tissue_interactions.index, 'normal_score'] += normal_only
    all_interactions.loc[tissue_interactions.index, 'tumor_score'] += tumor_only

n_tissues = len(tissues_cat) // 3

# Remove interactions that only have 0 values
keep = all_interactions[(all_interactions['both_score'] > 0) | (all_interactions['normal_score'] > 0) | (all_interactions['tumor_score'] > 0)].index
all_interactions = all_interactions.loc[keep]
all_interactions['diff'] = all_interactions['tumor_score'] - all_interactions['normal_score']

# Sort by difference between tumor and normal
all_interactions = all_interactions.sort_values('diff', ascending=False)
# Only save scores
all_interactions = all_interactions[['all_genes', 'n_genes', 'both_score', 'normal_score', 'tumor_score']]
all_interactions.to_csv('/home/lnemati/pathway_crosstalk/results/flow/bootstrap/bootstrap_scores.csv', index=True)

print("Done: aggregate_bootstrap.py")
