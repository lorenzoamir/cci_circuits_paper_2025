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
sub_tissue = filename.split("/")[-2]

print("Major tissue: {}".format(major_tissue))
print("Sub tissue: {}".format(sub_tissue))

# DEBUG
sys.exit()

# ----- Network ------

adj = WGCNA.adjacency # Pandas DataFrame

# TODO: read precomputed clusterings and subset to only genes that are in the network
# ORDER IS IMPORTANT

# Remove self-loops by filling diagonal with 0
#np.fill_diagonal(adj.values, 0)

# Create an igraph graph from the adjacency matrix
# Remove weak interactions that sum up to 10% of the total weight
# TODO: you can probably remove the 10% cutoff

# Flatten the matrix
flat_adj = adj.values
flat_adj = flat_adj[flat_adj > 0]
flat_adj = np.sort(flat_adj)

# Get total and cutoff values
total_weight = np.sum(flat_adj)
cutoff = 0.10 * total_weight

# Get weakest link to keep
cum_sum = np.cumsum(flat_adj)
num_to_remove = np.searchsorted(cum_sum, cutoff, side='right')
weakest_link = flat_adj[num_to_remove]

g = ig.Graph.Weighted_Adjacency(csr_matrix(np.where(adj < weakest_link, 0, adj)), mode=ig.ADJ_UNDIRECTED)

for clustering in all_clusterings:
    # Calculate network attributes
    modularity = g.modularity(clustering, weights=g.es['weight'])

    # TODO: save results, subtissue and major tissue should be retrievable

print("Done: modularity.py")
