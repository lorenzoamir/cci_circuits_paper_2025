import numpy as np
import sys
import scanpy as sc
import pandas as pd
import pingouin as pg
import os
import argparse

parser = argparse.ArgumentParser(description='Run WGCNA pipeline')

parser.add_argument('-i', '--input', type=str, help='path to input adata file')

args = parser.parse_args()

filename = args.input
print("Reading adata file: {}".format(filename))
adata = sc.read_h5ad(filename)

# name is the last part of the filename after the last /
name = filename.split("/")[-1].split(".")[0]
# output_path is the path to the name folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("name: ", name)
print("output_path: ", output_path)

# --------- WGCNA ------------

print("Starting WGCNA")

# make expression dataframe 
geneExp = pd.DataFrame(data=adata.X,  index=adata.obs_names, columns=adata.var_names)

# calculate the regular correlation matrix
corr_matrix = geneExp.corr()
outfile = os.path.join(output_path, 'correlation.csv.gz')
corr_matrix.to_csv(outfile, index=True, header=True, compression='gzip')

# calculate the partial correlation matrix using pingouin
partial_corr_matrix = pg.pcorr(geneExp)
outfile = os.path.join(output_path, 'partial_correlation.csv.gz')
partial_corr_matrix.to_csv(outfile, index=True, header=True, compression='gzip')

print("Done: partial_corr.py")
