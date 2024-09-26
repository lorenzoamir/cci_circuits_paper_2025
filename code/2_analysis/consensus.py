import PyWGCNA
import pandas as pd
import numpy as np
import gseapy as gp
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
    
    # Use one-hot encoding to get binary matrix of gene-module membership
    df = pd.get_dummies(wgcna.name + '_' + wgcna.datExpr.var.moduleLabels.astype(str))
    
    # Update genes set with the new genes
    genes.update(df.index)
    
    # Reindex the DataFrame to ensure all genes are present, filling missing values with 0
    combined_df = pd.concat([combined_df, df], axis=1).fillna(0)
    # Print the number of unique genes detected so far
    print(f"Number of unique genes detected: {len(genes)}")

# After the loop, combined_df will contain all genes with 0s for missing genes
df = combined_df
X = df.values

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.metrics import jaccard_score
import matplotlib.pyplot as plt

# Sample Sparse Binary Matrix (in a dense format for illustration)
# Each row is a binary vector

def jaccard_similarity_matrix(data):
    # Compute the intersection and union for Jaccard similarity
    intersection = np.dot(data, data.T)  # Number of 1s in both vectors
    union = data.sum(axis=1).reshape(-1, 1) + data.sum(axis=1) - intersection  # Total 1s in either vector

    # Compute Jaccard similarity
    similarity = intersection / union
    
    # Handle cases where union is zero (to avoid division by zero)
    similarity[np.isnan(similarity)] = 0  # Set NaN values to 0 (for all-zero rows)
    
    return similarity

def similarity_to_distance(similarity):
    return 1 - similarity  # Convert similarity to distance

def hierarchical_clustering(data):
    # Calculate Jaccard similarity
    print('Calculating similarity')
    jaccard_sim = jaccard_similarity_matrix(data)
    # Convert to distance matrix
    print('Calculating distance')
    distance_matrix = similarity_to_distance(jaccard_sim)
    
    # Perform hierarchical clustering using the 'ward' method
    # Note: You can change the method (e.g., 'average', 'single', 'complete')
    print('Calculating linkage')
    Z = linkage(distance_matrix, method='average')
   
    return Z

# Execute the hierarchical clustering
Z = hierarchical_clustering(df.values)

# Save the linkage matrix
linkage_matrix = pd.DataFrame(Z, columns=['cluster_1', 'cluster_2', 'distance', 'n_genes'])
linkage_matrix.to_csv(os.path.join(output, 'linkage_matrix.csv'), index=False)

# Save genes
genes_df = pd.DataFrame(list(genes), columns=['gene'])
genes_df.to_csv(os.path.join(output, 'genes.csv'), index=False)

print('Done: consensus.py ({})'.format(condition))
