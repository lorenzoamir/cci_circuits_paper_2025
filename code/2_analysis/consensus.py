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

print('List of WGCNA files:')
for path in file_list:
    print(path)

print('Reading all files')

combined_df = pd.DataFrame()

for path in file_list:
    print(path)

    # Read WGCNA object
    wgcna = PyWGCNA.readWGCNA(path)
    name = wgcna.name.replace('wgcna_', '')
    
    # Get KME module memeberships
    wgcna.CalculateSignedKME()
    df = wgcna.signedKME
    df.columns = [f'{name}_{x}' for x in df.columns]

    # Update genes set with the new genes
    genes.update(df.index)
    
    # Reindex the DataFrame to ensure all genes are present, use nans for missing genes
    combined_df = pd.concat([combined_df, df], axis=1)

    # Print the number of unique genes detected so far
    print(f"Number of unique genes detected: {len(genes)}")

# After the loop, combined_df will contain all genes with 0s for missing genes
df = combined_df.copy()
X = df.values

# Save the module membership matrix
df.to_csv(os.path.join(output, 'module_membership.csv'))

# Create correlation matrix
correlation_matrix = df.T.corr()
correlation_matrix.to_csv(os.path.join(output, 'correlation_matrix.csv'))

distance_matrix = 1 - correlation_matrix

# Execute the hierarchical clustering
Z = linkage(distance_matrix, method='average')

# Save the linkage matrix
linkage_matrix = pd.DataFrame(Z, columns=['cluster_1', 'cluster_2', 'distance', 'n_genes'])
linkage_matrix.to_csv(os.path.join(output, 'linkage_matrix.csv'), index=False)

# Generate all possible clusterings
genes = pd.DataFrame(list(genes), columns=['gene']).set_index('gene')

for i in range(2, 100):
    genes[f'{i}_clusters'] = fcluster(Z, i, criterion='maxclust')

# Save genes
genes.to_csv(os.path.join(output, 'all_clusterings.csv'))

print('Done: consensus.py ({})'.format(condition))
