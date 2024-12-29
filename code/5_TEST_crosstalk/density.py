import pandas as pd
import os
import numpy as np
import gseapy as gp
import argparse
import igraph as ig
from igraph import Graph
from itertools import combinations
from multiprocessing import Pool
import ast
import sys
import re

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--motifs', type=str, help='Path to motifs.csv file')
parser.add_argument('--network', type=str, help='Path to co-expression network file')

args = parser.parse_args()

motifs_path = args.motifs
network_path = args.network

network = pd.read_csv(
    network_path,
    usecols=['complex1', 'complex2', 'adj']
)

# Add a mirrored version of the network to make it symmetric
network_mirrored = network.copy()
network_mirrored.columns = ['complex2', 'complex1', 'adj']
network = pd.concat([network, network_mirrored], ignore_index=True)

adj = network.pivot(index='complex1', columns='complex2', values='adj').fillna(0)

def get_density(interaction, adj):
    # Get the density (sum of weights / N possible edges)
    # for the nodes in the interaction

    # Get the nodes in the interaction
    nodes = re.split('[+&]', interaction)
    n_nodes = len(nodes)
    n_edges = n_nodes * (n_nodes - 1) / 2

    # Subset adj to the nodes in the interaction
    adj_sub = adj.loc[nodes, nodes]
    # Get upper triangle of adj_sub
    adj_sub = np.triu(adj_sub, k=1)
    # Get sum of weights
    sum_weights = adj_sub.sum()
    # Get density
    density = sum_weights / n_edges
    
    return density

motifs = pd.read_csv(motifs_path)

# Use multiprocessing to get the density for each motif
ncpus = int(os.getenv('NCPUS', default=1))
pool = Pool(ncpus) 
densities = pool.starmap(get_density, [(motif, adj) for motif in motifs['interaction']])
pool.close()

# Add the densities to the motifs dataframe
motifs['density'] = densities

# Sort based on motif type and then density
motifs = motifs.sort_values(by=['motif', 'density'], ascending=[True, False])
motifs.to_csv(motifs_path, index=False)

print('Done: density.py')
