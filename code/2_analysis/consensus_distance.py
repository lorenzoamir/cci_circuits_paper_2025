import PyWGCNA
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import cosine_distances
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate linkage for consensus modules')

parser.add_argument('--file-list', type=str, help='List of directories to compare, delimited by commas')
parser.add_argument('--condition', type=str, help='tumor or normal')

args = parser.parse_args()

file_list = args.file_list
condition = args.condition

# Only keep directories that contain the condition
file_list = file_list.split(',')
file_list = [x.strip() for x in file_list]
file_list = [d for d in file_list if '/{}/'.format(condition) in d]

output = '/home/lnemati/pathway_crosstalk/results/consensus_modules/{}'.format(condition)
# Create output directory if it does not exist
os.makedirs(output, exist_ok=True)

genes = set()

# Get all genes by looking at all wgcna files
for path in file_list:
    wgcna = PyWGCNA.readWGCNA(path)
    genes.update(wgcna.datExpr.var.index)

genes = list(genes)
triu_indices = np.triu_indices(len(genes), k=1)
print('Found {} unique genes'.format(len(genes)))

flattened = []

print('Reading all files')

def get_all_toms(paths, genes):
    toms = []
    triu_indices = np.triu_indices(len(genes), k=1)
    for path in paths:
        wgcna = PyWGCNA.readWGCNA(path)
        tom = wgcna.TOM
        tom = tom.reindex(index=genes, columns=genes, fill_value=0)
        toms.append(tom.values[triu_indices])
    toms = np.vstack(toms)
    return toms

def get_all_adjs(paths, genes):
    adjs = []
    triu_indices = np.triu_indices(len(genes), k=1)
    for path in paths:
        wgcna = PyWGCNA.readWGCNA(path)
        adj = wgcna.adjacency
        # Normalize adj so that max degree is 1
        degree = adj.sum(axis=1)
        adj = adj / max(degree)
        adj = adj.reindex(index=genes, columns=genes, fill_value=0)
        adjs.append(adj.values[triu_indices])
    adjs = np.vstack(adjs)
    return adjs

def get_all_corrs(paths, genes):
    corrs = []
    triu_indices = np.triu_indices(len(genes), k=1)
    for path in paths:
        wgcna = PyWGCNA.readWGCNA(path)
        # calculate correlation matrix
        corr = pd.DataFrame(np.corrcoef(wgcna.datExpr.X.T), index=wgcna.datExpr.var.index, columns=wgcna.datExpr.var.index)
        # reindex to include all genes
        corr = corr.reindex(index=genes, columns=genes, fill_value=0)
        corrs.append(corr.values[triu_indices])
    corrs = np.vstack(corrs)
    return corrs

for metric, function in [('tom', get_all_toms), ('adj', get_all_adjs), ('corr', get_all_corrs)]:
    print('Calculating {}'.format(metric))
    stacked = function(file_list, genes)
    # Calculate median across tissues (one value per gene pair)
    print('Calculating median')
    flat = np.median(stacked, axis=0)
    # Reshape to square matrix
    median = np.zeros((len(genes), len(genes)))
    median[triu_indices] = flat
    median = median + median.T
    median = pd.DataFrame(median, index=genes, columns=genes)
    median.to_csv(os.path.join(output, 'median_{}.csv'.format(metric)))

print('Done: consensus.py ({})'.format(condition))

## !!! This is the old code before converting to functions
#for path in file_list:
#    print(path)
#
#    # Read WGCNA object
#    wgcna = PyWGCNA.readWGCNA(path)
#    tom = wgcna.TOM
#    
#    # Expand tom (both rows and cols) to include all genes, use 0 for missing genes
#    tom = tom.reindex(index=genes, columns=genes, fill_value=0)
#
#    # Only use upper triangle to save memory
#    flattened.append(tom.values[triu_indices])

## Calculate median across tissues (one value per gene pair)
#print('Calculating median')
#flat = np.median(stacked, axis=0)
## Reshape to square matrix
#median = np.zeros((len(genes), len(genes)))
#median[triu_indices] = flat
#median = median + median.T
#median = pd.DataFrame(median, index=genes, columns=genes)
#clustering_pipeline(median, os.path.join(output, 'median'))
#
## Do the same with 25% percentile
#print('Calculating 25th percentile')
#flat = np.percentile(stacked, 25, axis=0)
## Reshape to square matrix
#perc25 = np.zeros((len(genes), len(genes)))
#perc25[triu_indices] = flat
#perc25 = perc25 + perc25.T
#perc25 = pd.DataFrame(perc25, index=genes, columns=genes)
#clustering_pipeline(perc25, os.path.join(output, 'perc25'))
#

#def clustering_pipeline(matrix, path, max_clusters=250):
#    os.makedirs(path, exist_ok=True)
#    # Save the matrix
#    matrix.to_csv(os.path.join(path, 'matrix.csv'))
#    genes = matrix.index
#    # Calculate dissimilarity matrix
#    dissimilarity = 1 - matrix
#    dissimilarity.to_csv(os.path.join(path, 'dissimilarity.csv'))
#    dissimilarity = squareform(dissimilarity)
#    # Generate linkage
#    Z = linkage(dissimilarity, method='average')
#    Z = pd.DataFrame(Z, columns=['cluster_1', 'cluster_2', 'distance', 'n_genes'])
#    Z.to_csv(os.path.join(path, 'linkage.csv'), index=False)
#
#    return
#
