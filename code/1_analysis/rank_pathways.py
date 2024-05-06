import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Rank high level pathways based on links between low level pathways')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

pw_net= pd.read_csv('/home/lnemati/pathway_crosstalk/results/pathways_network/pathways_network_filtered.csv')
roots = pd.read_csv('/home/lnemati/resources/reactome/ReactomeRootPathways.csv')

print('Pathway network:')
print(pw_net.head())

# Read .gmt to get association between pathway id and name
print('Reading .gmt file')
pw_dict = {}

with open('/home/lnemati/resources/reactome/ReactomePathways.gmt') as f:
    lines = f.readlines()
    for line in lines:
        line = line.strip().split('\t')
        pw_dict[line[1]] = line[0]

print('ID to name dictionary:')
for key in list(pw_dict.keys())[:3]:
    print(key, pw_dict[key])

print('Initializing adjacency matrix')
# Make adjacency matrix with root pathways as rows and columns
adj = pd.DataFrame(0, index=roots['root'].unique(), columns=roots['root'].unique())

print('Shape: ', adj.shape)
print('Adjacency matrix:')
print(adj.head())

print('Filling adjacency matrix')

# Fill adjacency matrix
for root_id1 in adj.index:
    for root_id2 in adj.columns:
        # Get all children
        children1 = roots[roots['root'] == root_id1]['pathway']
        children1 = [pw_dict[child] for child in children1 if child in pw_dict.keys()]
        children2 = roots[roots['root'] == root_id2]['pathway']
        children2 = [pw_dict[child] for child in children2 if child in pw_dict.keys()]

        # Get all edges between children
        edges = pw_net[(pw_net['pathway1'].isin(children1)) & (pw_net['pathway2'].isin(children2))]
        adj.loc[root_id1, root_id2] = len(edges)

# Make symmetric
print('Making symmetric')
adj = adj + adj.T

# Convert root ids to names
print('Converting root ids to names')
adj.columns = [pw_dict[col] for col in adj.columns]
adj.index = [pw_dict[ind] for ind in adj.index]

# Subset to only include pathways with at least one edge
print('Subsetting to only include pathways with at least one edge')
adj = adj.loc[adj.sum(axis=1) > 0, adj.sum(axis=0) > 0]



# Save adjacency matrix
print('Saving adjacency matrix')
adj.to_csv('/home/lnemati/pathway_crosstalk/results/pathways_network/pathways_high_level_adjacency.csv')

print('Done: rank_pathways.py')
