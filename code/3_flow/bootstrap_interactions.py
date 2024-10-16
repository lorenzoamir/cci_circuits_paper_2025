import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import ast
import os
import argparse

parser = argparse.ArgumentParser(description='Find interactions in WGCNA results')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')

args = parser.parse_args()

filename = args.input
print('Reading WGCNA file:', filename)
WGCNA = PyWGCNA.readWGCNA(filename)

# category is the last part of the filename after the last /
category = filename.split("/")[-1].split(".")[0]
# output_path is the path to the category folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("category", category)
print("output_path", output_path)
print("filename", filename)

# --------- interactions ---------

genes = WGCNA.datExpr.var.index
print("Number of genes in WGCNA data: ", len(genes))
print("Showing the first 10 genes")
print(genes[:10])

bootstrap_path = '/home/lnemati/pathway_crosstalk/results/flow/bootstrap'

# Read bootstap_all.csv
print("Reading bootstrap_all.csv")
interactions = pd.read_csv(os.path.join(bootstrap_path, "bootstrap_all.csv"), index_col=0)

# The all_genes column contains a list of genes and should not be interpreted as a string
# Remove first and last character ('[' and ']') and split by ', ' to get a list of genes
interactions["all_genes"] = interactions["all_genes"].apply(lambda x: ast.literal_eval(x))
print("Showing the first 10 interactions")
print(interactions.head(10))

# Only keep interactions where all genes are in the WGCNA data
print(f'Number of interactions before filtering: {len(interactions)}')
interactions = interactions[interactions.all_genes.apply(lambda x: all([gene in genes for gene in x]))]
print("Keeping interactions with all genes present: ", len(interactions))

# Remove interactions with only one gene
interactions = interactions[interactions.all_genes.apply(lambda x: len(set(x)) > 1)]
print("Keeping interactions with more than one gene: ", len(interactions))

result = interactions.copy()
result["same_module"] = False

for i, row in interactions.iterrows():
    all_genes = row["all_genes"]

    if WGCNA.datExpr.var.loc[all_genes, "moduleLabels"].nunique() == 1:
        result.loc[i, "same_module"] = True
    
print("Showing the first 10 interactions")
print(result.head(10))

# Save the interactions in the same path as the WGCNA object
print("Saving interactions to {}".format(os.path.join(output_path, "bootstrap_interactions.csv")))
result.to_csv(os.path.join(output_path, "interactions_bootstrap.csv"), index=True)

print("Done: bootstrap_interactions.py")
