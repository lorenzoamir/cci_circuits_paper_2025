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


# Get strength of connections between hubs and pathways
results = pd.DataFrame(index=root_pws, columns=['connection_strength'])

for pw in root_pws:
    print(pw)
    genes = list(set(pw_genes[pw]).intersection(adj.columns))
    tot = np.sum(adj.loc[hubs, genes].values)
    results.loc[pw, 'connection_strength'] = tot

# save degree_df as csv
results.to_csv(os.path.join(output_path, "hubs_connectivities.csv"))

print("Done: hubs_connectivities.py")
