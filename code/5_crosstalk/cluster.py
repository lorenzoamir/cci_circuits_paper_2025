import pandas as pd
import os
import numpy as np
import gseapy as gp
import argparse
import igraph as ig
from igraph import Graph
import leidenalg as la
import ast
import sys

# Read unfiltered network
network = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/network_unfiltered.csv.gz')

# Extract edges (node pairs) and weights
edges = list(zip(network['complex1'], network['complex2']))

# Renormalize AUROC to [-1, 1], 1 = tumor more co-expressed, -1 = normal more co-expressed
weights = (2*network['auroc'] - 1).tolist()

# Create the undirected graph
G = Graph.TupleList(edges, directed=False, edge_attrs=['weight'])
G.es['weight'] = weights

# Also get adjacency matrix
matrix = pd.DataFrame(G.get_adjacency(attribute='weight'), index=G.vs['name'], columns=G.vs['name'])

# Create positive and negative layers
G_pos = G.subgraph_edges(G.es.select(weight_gt = 0), delete_vertices=False);
G_neg = G.subgraph_edges(G.es.select(weight_lt = 0), delete_vertices=False);
# Reverse sign of negative layer
G_neg.es['weight'] = [-w for w in G_neg.es['weight']];

# Run Leiden algorithm
partition = la.find_partition_multiplex(
    graphs = [G_pos, G_neg],
    partition_type = la.ModularityVertexPartition,
    layer_weights = [1, -1],
    #max_comm_size = 500,
    n_iterations=-1,
    weights='weight',
    seed=42,
)

clusters = pd.Series(partition[0], index=G.vs['name'])

print('Number of clusters with more than 1 complex:')
print(clusters.value_counts()[clusters.value_counts() > 1])

print('Number of complexes in each cluster:')
clusters.value_counts()

# Save clusters
outdir = '/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/clusters'
os.makedirs(outdir, exist_ok=True)





print('Done: cluster.py')
