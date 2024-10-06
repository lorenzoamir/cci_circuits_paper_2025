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
for path in file_list:
    print(path)

    # Read WGCNA object
    wgcna = PyWGCNA.readWGCNA(path)
    corr = np.corrcoef(wgcna.datExpr.X.T)

    # Expand corr (both rows and cols) to include all genes, use 0 for missing genes
    corr = pd.DataFrame(corr, index=wgcna.datExpr.var.index, columns=wgcna.datExpr.var.index)
    corr = corr.reindex(index=genes, columns=genes, fill_value=0)

    # Only use upper triangle to save memory
    flattened.append(corr.values[triu_indices])

print()
# Stack all vectors into a 2D numpy array
print('Stacking matrices')
stacked = np.vstack(flattened)

# Calculate median across tissues (one value per gene pair)
print('Calculating median')
flat = np.median(stacked, axis=0)
# Reshape to square matrix
median = np.zeros((len(genes), len(genes)))
median[triu_indices] = flat
median = median + median.T
median = pd.DataFrame(median, index=genes, columns=genes)
# Save the matrix
median.to_csv(os.path.join(output, 'median', 'correlation.csv')) 

print('Done: correlation_consensus.py ({})'.format(condition))
