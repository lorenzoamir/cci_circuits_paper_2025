import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
from scipy.stats import spearmanr
import os
import argparse

parser = argparse.ArgumentParser(description='Network analysis of WGCNA network')

parser.add_argument('-d', '--directory', type=str, help='path to wgcna dir, serves as input path for tom and output path')
parser.add_argument('-i', '--input', type=str, help='path to input coevoulution matrix')

args = parser.parse_args()

directory = args.directory
coevolution = args.input

# ------ TOM & adjacency ------

print("reading TOM matrix from: ", os.path.join(directory, "tom.csv.gz"))
tom = pd.read_csv(os.path.join(directory, "tom.csv.gz"), index_col=0)

print("Reading coevolution matrix from: ", coevolution)
coevolution = pd.read_csv(coevolution, index_col=0)

genes = tom.index.intersection(coevolution.index)

print("Genes in TOM: ", tom.shape[0])
print("Genes in coevolution: ", coevolution.shape[0])
print(f"Genes in common: {genes.shape[0]}, ({genes.shape[0] / tom.shape[0] * 100:.2f}%)")

# Raise warning if less than 75% of genes are in common
if genes.shape[0] / tom.shape[0] < 0.75:
    print("Less than 75% of genes are in common between TOM and coevolution matrix")
    print("Less than 75% of genes are in common between TOM and coevolution matrix", file=sys.stderr)

# Subset and reorder matrices
tom = tom.loc[genes, genes]
coevolution = coevolution.loc[genes, genes]

# Correlation
correlation = pd.DataFrame(index=genes, columns=['spearman_r', 'pvalue'])

for gene in tom.index:
    # exclude diagonal
    tom_values = tom.loc[gene, tom.columns != gene].values
    coevolution_values = coevolution.loc[gene, coevolution.columns != gene].values
    correlation.loc[gene, 'spearman_r'], correlation.loc[gene, 'pvalue'] = spearmanr(tom_values, coevolution_values)

# Save correlation as csv
correlation.to_csv(os.path.join(directory, "tom_coevolution_correlation.csv"))

average_correlation = correlation['spearman_r'].mean()

print("Writing general info")
# Append, don't overwrite
with open(os.path.join(directory, "general_info.txt"), 'a') as f:
    f.write("average_tom_coev_correlation: {}\n".format(average_correlation))

print("Done: coevolution.py")
