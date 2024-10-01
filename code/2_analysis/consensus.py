import PyWGCNA
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
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
print('Found {} unique genes'.format(len(genes)))

matrices = []

print('Reading all files')
for path in file_list:
    print(path)

    # Read WGCNA object
    wgcna = PyWGCNA.readWGCNA(path)
    tom = wgcna.TOM
    
    # Expand tom (both rows and cols) to include all genes, use 0 for missing genes
    tom = tom.reindex(index=genes, columns=genes, fill_value=0)
   
    # Add the matrix to the list
    matrices.append(tom.values)

print()
# Stack all matrices into a 3D numpy array
print('Stacking matrices')
stacked = np.dstack(matrices)

def clustering_pipeline(matrix, path, max_clusters=250):
    os.makedirs(path, exist_ok=True)
    # Save the matrix
    matrix.to_csv(os.path.join(path, 'matrix.csv'))
    genes = matrix.index
    # Calculate dissimilarity matrix
    dissimilarity = 1 - matrix
    dissimilarity.to_csv(os.path.join(path, 'disimilarity.csv'))
    # Generate linkage
    Z = linkage(dissimilarity, method='average')
    Z = pd.DataFrame(Z, columns=['cluster_1', 'cluster_2', 'distance', 'n_genes'])
    Z.to_csv(os.path.join(path, 'linkage.csv'), index=False)
    # Generate possible clusterings
    all_cluterings = pd.DataFrame(list(genes), columns=['gene']).set_index('gene')
    for i in range(2, max_clusters):
        all_cluterings[f'{i}_clusters'] = fcluster(Z, i, criterion='maxclust')
    all_cluterings.to_csv(os.path.join(path, 'all_clusterings.csv'))

    return

# Calculate median
print('Calculating median')
median = np.median(stacked, axis=2)
median = pd.DataFrame(median, index=genes, columns=genes)
clustering_pipeline(median, os.path.join(output, 'median'))

# Do the same with 25% percentile
perc25 = np.quantile(stacked, 0.25, axis=2)
perc25 = pd.DataFrame(perc25, index=genes, columns=genes)
clustering_pipeline(perc25, os.path.join(output, 'perc25'))

print('Done: consensus.py ({})'.format(condition))
