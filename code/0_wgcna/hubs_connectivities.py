import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Network analysis of WGCNA network')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')

args = parser.parse_args()

# Read input file
filename = args.input
output_path = os.path.dirname(filename)
WGCNA = PyWGCNA.readWGCNA(filename)
adj = WGCNA.adjacency
np.fill_diagonal(adj.values, 1)

# Get the top 5% hubs
hubs = list(adj.sum().sort_values(ascending=False).iloc[:round(0.05 * len(adj))].index)

# Read pathways
pw_genes   = {}
name_to_id = {}
id_to_name = {}

with open('/home/lnemati/resources/reactome/ReactomePathways.gmt', 'r') as f:
    for line in f.readlines():
        line = line.split('\t')
        name_to_id[line[0]] = line[1]
        id_to_name[line[1]] = line[0]
        pw_genes[line[0]] = line[2:]

root_pws = pd.read_csv('/home/lnemati/resources/reactome/ReactomeRootPathways.csv')
root_pws = root_pws['root'].map(id_to_name).unique()
all_pws  = pw_genes.keys()

print('Getting total strength of connections between hubs and pathways')
total = adj.loc[hubs, :].sum().sum()

# Normalize adjacency matrix so that the sum of 
# all connections of hubs to pathways is 1
print('Normalizing adjacency matrix')
adj = adj / total

# Get strength of connections between hubs and pathways
print('Getting strength of connections between hubs and each pathway')
results = pd.DataFrame(index=root_pws, columns=['connection_strength'])

# Calculate connection strength with root pathways
for pw in root_pws:
    print(pw)
    genes = list(set(pw_genes[pw]).intersection(adj.columns))
    tot = np.sum(adj.loc[hubs, genes].values)
    results.loc[pw, 'connection_strength'] = tot

# save degree_df as csv
results.to_csv(os.path.join(output_path, "hubs_connectivities_root.csv"))

# Calculate connection strength with all pathways
for pw in all_pws:
    print(pw)
    genes = list(set(pw_genes[pw]).intersection(adj.columns))
    tot = np.sum(adj.loc[hubs, genes].values)
    results.loc[pw, 'connection_strength'] = tot

# save degree_df as csv
results.to_csv(os.path.join(output_path, "hubs_connectivities_all.csv"))

# Also do the same with MSigDB Hallmark genesets
hallmarks = '/home/lnemati/resources/msigdb_hallmarks/h.all.v2024.1.Hs.symbols.gmt'
hallmark_genes = {}

with open(hallmarks, 'r') as f:
    for line in f.readlines():
        line = line.split('\t')
        hallmark_genes[line[0]] = line[2:]

# Init new results dataframe
results = pd.DataFrame(index=list(hallmark_genes.keys()), columns=['connection_strength'])

# Calculate connection strength with MSigDB Hallmark pathways
for pw in hallmark_genes.keys():
    print(pw)
    genes = list(set(hallmark_genes[pw]).intersection(adj.columns))
    tot = np.sum(adj.loc[hubs, genes].values)
    results.loc[pw, 'connection_strength'] = tot

# save degree_df as csv
results.to_csv(os.path.join(output_path, "hubs_connectivities_hallmarks.csv"))

print("Done: hubs_connectivities.py")
