import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from dynamicTreeCut import cutreeHybrid
import sys
import os
import argparse

parser = argparse.ArgumentParser(description='Find consensus modules from consensus tom matrix')

parser.add_argument('-i', '--inputfile', type=str, help='path to tom matrix', required=True)

args = parser.parse_args()

inputfile = args.inputfile
inputdir = os.path.dirname(inputfile)

#link = pd.read_csv(os.path.join(inputdir, 'linkage.csv'), index_col=0).values
distances = 1 - pd.read_csv(os.path.join(inputfile), index_col=0)
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
