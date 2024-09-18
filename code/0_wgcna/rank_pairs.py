import numpy as np
import sys
import PyWGCNA
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Network analysis of WGCNA network')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')

args = parser.parse_args()

filename = args.input
WGCNA = PyWGCNA.readWGCNA(filename)

# category is the last part of the filename after the last /
category = filename.split("/")[-1].split(".")[0]
# output_path is the path to the category folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("category", category)
print("output_path", output_path)
print("filename", filename)

# ------ Rank adjacency ------

print("Getting adjacency matrix")

def rank_matrix(matrix):
    original_shape = matrix.shape

    # Only keep upper triangle, set diagonal and lower triangle to 0
    matrix = np.triu(matrix , 1)

    n_elements = matrix[matrix > 0].size
    print(f'n_elements: {n_elements}')

    # Flatten the matrix into a 1D array
    flat_matrix = matrix.flatten()

    # Get the ranks (argsort returns the indices that would sort the array)
    ranks = np.argsort(np.argsort(flat_matrix))

    print(f'max rank: {ranks.max()}')

    ranks = ranks.reshape(original_shape)
    print(ranks)

    # Use n_elements to normalize the ranks so that all elements outside the upper triangle are 0
    threshold = np.max(ranks) - n_elements
    print(f'threshold: {threshold}')

    ranks = np.where(ranks > threshold, ranks - threshold, 0)
    print(ranks)
    print()

    print(f'new max: {ranks.max()}')
    print(ranks)
    assert (ranks == np.triu(ranks , 1)).all()

    print(f'non 0: {np.sum(ranks > 0)}')
    assert (ranks == np.triu(ranks , 1)).all()

    # Reshape the ranks back to the original matrix shape
    ranks = ranks.reshape(original_shape)
    print(f'new max: {ranks.max()}')
    assert (ranks == np.triu(ranks , 1)).all()

    # Normalizing
    print('Normalizing')
    ranks = ranks / n_elements
    print(f'new max: {ranks.max()}')
    assert (ranks.max() == 1)
    # Make symmetric

    print('Making symmetric')
    ranks = ranks + ranks.T
    print(f'new max: {ranks.max()}')
    assert (ranks.max() == 1)

    return ranks

print("Getting adjacency matrix")
adj = WGCNA.adjacency
index = adj.index

# Rank adjacency matrix
print("Ranking adjacency matrix")
# Back to dataframe
ranked_adj = pd.DataFrame(rank_matrix(adj), index=index, columns=index)
ranked_adj.to_csv(os.path.join(output_path, "ranked_adjacency.csv.gz"), compression="gzip")

print("Done: network.py")
