import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from dynamicTreeCut import cutreeHybrid
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Find consensus modules from dissimilarity matrix')

parser.add_argument('-i', '--inputdir', type=str, help='path to directory containing dissimilarity.csv', required=True)

args = parser.parse_args()

inputdir = args.inputdir

#link = pd.read_csv(os.path.join(inputdir, 'linkage.csv'), index_col=0).values
distances = pd.read_csv(os.path.join(inputdir, 'dissimilarity.csv'), index_col=0)
# Fill diagonal with 0s
np.fill_diagonal(distances.values, 0)

print('Distance matrix shape:', distances.shape)
print(distances.head())

genes = distances.index

distances = squareform(distances.values)
link = linkage(distances, method='average')

clusters = cutreeHybrid(link, distances, minClusterSize=20)
module_genes = pd.Series(clusters['labels'], index=genes, name='module')

# Save to csv, the column is called 'module'
module_genes.to_csv(os.path.join(inputdir, 'consensus_modules.csv'))

print("Done: cluster.py")
