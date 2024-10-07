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
parser.add_argument('--condition', type=str, help='condition', required=True)
parser.add_argument('--quantile', type=str, help='quantile', required=True)

args = parser.parse_args()

filename = args.input
WGCNA = PyWGCNA.readWGCNA(filename)
condition = args.condition
quantile = args.quantile

output_path = '/home/lnemati/pathway_crosstalk/results/consensus_modules'

# category is the last part of the filename after the last /
major_tissue = filename.split("/")[-4]
#condition = filename.split("/")[-3]
sub_tissue = filename.split("/")[-2]

# Check if condition is in the filename, else exit
if filename.split("/")[-3] != condition:
    print("Condition does not match the filename")
    print("Filename: {}".format(filename))
    print("Condition: {}".format(condition))
    print("Done: modularity.py")
    sys.exit(0)

print("Condition: {}".format(condition))
print("Quantile: {}".format(quantile))
print("Major tissue: {}".format(major_tissue))
print("Sub tissue: {}".format(sub_tissue))

# ----- Network ------

adj = WGCNA.adjacency # Pandas DataFrame
adj = (adj + adj.T) / 2

# Read precomputed clusterings
print("Reading consensus modules")
modules_path = f"/home/lnemati/pathway_crosstalk/results/consensus_modules/{condition}/{quantile}/consensus_modules.csv"

modules = pd.read_csv(modules_path, index_col=0)
modules = modules['module']

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

out_tissue_path = os.path.join(output_path, condition, quantile, 'tissues', major_tissue, sub_tissue)

# Calculate modularity
clustering = modules.loc[adj.index].values
modularity = g.modularity(clustering, weights=g.es['weight'])

# Create modularity file in the tissue folder
os.makedirs(out_tissue_path, exist_ok=True)
modularity_path = os.path.join(out_tissue_path, 'modularity.txt')
with open(modularity_path, 'w') as f:
    f.write(str(modularity))

print("Done: modularity.py")
