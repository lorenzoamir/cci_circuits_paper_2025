import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import powerlaw
import os
import igraph as ig
from scipy.sparse import csr_matrix
import argparse

parser = argparse.ArgumentParser(description='Network analysis of WGCNA network')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')

args = parser.parse_args()

filename = args.input
WGCNA = PyWGCNA.readWGCNA(filename)

output_path = '/home/lnemati/pathway_crosstalk/results/consensus_modules'

# category is the last part of the filename after the last /
major_tissue = filename.split("/")[-4]
condition = filename.split("/")[-3]
sub_tissue = filename.split("/")[-2]

print("Major tissue: {}".format(major_tissue))
print("Condition: {}".format(condition))
print("Sub tissue: {}".format(sub_tissue))

# ----- Network ------

adj = WGCNA.adjacency # Pandas DataFrame
adj = (adj + adj.T) / 2

# Read precomputed clusterings
print("Reading precomputed clusterings")

quantiles = ['perc25', 'median']

for quantile in quantiles:

    if quantile == 'perc25':
        # DEBUG
        print("Skipping perc25")
        continue

    clustering_path = os.path.join(output_path, condition, quantile, 'all_clusterings.csv')
    all_clusterings = pd.read_csv(clustering_path, index_col=0)

    # Remove self-loops by filling diagonal with 0
    #np.fill_diagonal(adj.values, 0)

    # Create an igraph graph from the adjacency matrix
    # Remove weak interactions that sum up to 10% of the total weight

    # Flatten the matrix
    #flat_adj = adj.values
    #flat_adj = flat_adj[flat_adj > 0]
    #flat_adj = np.sort(flat_adj)
    #
    ## Get total and cutoff values
    #total_weight = np.sum(flat_adj)
    #cutoff = 0.10 * total_weight
    #
    ## Get weakest link to keep
    #cum_sum = np.cumsum(flat_adj)
    #num_to_remove = np.searchsorted(cum_sum, cutoff, side='right')
    #weakest_link = flat_adj[num_to_remove]
    #
    #g = ig.Graph.Weighted_Adjacency(csr_matrix(np.where(adj < weakest_link, 0, adj)), mode=ig.ADJ_UNDIRECTED)

    g = ig.Graph.Weighted_Adjacency(adj.values, mode=ig.ADJ_UNDIRECTED)
    result = pd.DataFrame()

    out_tissue_path = os.path.join(output_path, condition, quantile, 'tissues', major_tissue, sub_tissue)

    for column in all_clusterings.columns:

        # Calculate modularity
        clustering = all_clusterings.loc[adj.index, column].values
        modularity = g.modularity(clustering, weights=g.es['weight'])

        # Use current column as index of result and add modularity
        result.at[column, 'modularity'] = modularity

        # Create dir and save the clustering
        os.makedirs(out_tissue_path, exist_ok=True)
        result.to_csv(os.path.join(out_tissue_path, 'modularity.csv'))

print("Done: modularity.py")
