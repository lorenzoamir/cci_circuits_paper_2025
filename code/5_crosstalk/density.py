import pandas as pd
import os
import numpy as np
import argparse
from multiprocessing import Pool
import ast
import sys
import re

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--motifs', type=str, help='Path to motifs.csv file')
args = parser.parse_args()

motifs_path = args.motifs
print('motifs_path:', motifs_path)

network_dir = '/home/lnemati/pathway_crosstalk/data/networks'

# Inside network_dir are normal and tumor directories, save their realpaths
normal_nets_paths = [os.path.join(network_dir, 'normal', f) for f in os.listdir(os.path.join(network_dir, 'normal'))]
tumor_nets_paths = [os.path.join(network_dir, 'tumor', f) for f in os.listdir(os.path.join(network_dir, 'tumor'))]
print('Normal networks:', normal_nets_paths)
print('Tumor networks:', tumor_nets_paths)

def network_to_adj(network_path):
    network = pd.read_csv(
        network_path,
        usecols=['complex1', 'complex2', 'adj']
    )

    # Add a mirrored version of the network to make it symmetric
    network_mirrored = network.copy()
    network_mirrored.columns = ['complex2', 'complex1', 'adj']
    network = pd.concat([network, network_mirrored], ignore_index=True)

    adj = network.pivot(index='complex1', columns='complex2', values='adj').fillna(0)
    return adj

def get_density(interaction, adj):
    # Get the density (sum of weights / N possible edges)
    # for the nodes in the interaction

    # Get the nodes in the interaction
    nodes = re.split('[+&]', interaction)
    n_nodes = len(nodes)

    # Get the number of possible edges before removing missing nodes
    n_edges = n_nodes * (n_nodes - 1) / 2

    # Subset adj to the nodes in the interaction
    nodes = adj.index.intersection(nodes)
    
    # If there are no nodes, return 0
    if len(nodes) == 0:
        return 0

    adj_sub = adj.loc[nodes, nodes]
    # Get upper triangle of adj_sub
    adj_sub = np.triu(adj_sub, k=1)
    # Get sum of weights
    sum_weights = adj_sub.sum()
    # Get density
    density = sum_weights / n_edges
    
    return density

motifs = pd.read_csv(motifs_path)

# Rename columns, Type -> motif, Interaction -> interaction
motifs = motifs.rename(columns={'Type': 'motif', 'Interaction': 'interaction'})

# Get condition (normal, tumor or both) from the directory name
condition = os.path.basename(os.path.dirname(motifs_path))
print('Condition:', condition)

ncpus = int(os.getenv('NCPUS', default=1))
print('ncpus:', ncpus)

# Get the network path based on the condition
network_paths = normal_nets_paths + tumor_nets_paths

print(f'Using {len(network_paths)} networks')

all_names = []
for network_path in network_paths:
    print('Network:', network_path)
    # Get tissue from filename
    tissue = os.path.basename(network_path)
    # Remove extension
    tissue = os.path.splitext(tissue)[0]
    # Set name to condition_tissue
    if network_path in normal_nets_paths:
        name = 'normal_' + tissue
    elif network_path in tumor_nets_paths:
        name = 'tumor_' + tissue
    all_names.append(name)

    # Read the network
    adj = network_to_adj(network_path)

    # Get the density for each motif
    pool = Pool(ncpus)
    densities = pool.starmap(get_density, [(motif, adj) for motif in motifs['interaction']])
    pool.close()

    # Save the densities to the motifs dataframe
    motifs[name] = densities

# Average the densities across all networks
motifs['avg_tumor'] = motifs[[name for name in all_names if 'tumor_' in name]].mean(axis=1)
motifs['avg_normal'] = motifs[[name for name in all_names if 'normal_' in name]].mean(axis=1)
motifs['avg_all'] = (motifs['avg_tumor'] + motifs['avg_normal']) / 2

# Reorder the columns: Type, Interaction, average_density, then all the network densities in alphabetical order
# make all columns lowercase, remove spaces
motifs.columns = motifs.columns.str.lower().str.replace(' ', '_')
all_names.sort()
motifs = motifs[['motif', 'interaction', 'avg_all', 'avg_normal', 'avg_tumor'] + all_names]

# Sort based on motif type and then density
if condition == 'both':
    sort_col = 'avg_all'
elif condition == 'normal':
    sort_col = 'avg_normal'
elif condition == 'tumor':
    sort_col = 'avg_tumor'

motifs = motifs.sort_values(by=['motif', sort_col], ascending=[True, False])

out_dir = motifs_path.replace('motifs.csv', '')
motifs.to_csv(os.path.join(out_dir, 'motifs_densities.csv'), index=False)

print('Done: density.py')
