import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Find interactions in WGCNA results')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')
parser.add_argument('-r', '--lr_resource', type=str, help='path to ligand-receptor resource')

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

print("Loading LR resource: {}".format(args.lr_resource))
interactions = pd.read_csv(args.lr_resource)

# Extract genes for complex A and B for each row
complex_a_genes = interactions.apply(lambda row: [row[f'interactor{i}'] for i in range(1, row['num_interactors_a'] + 1) if pd.notna(row[f'interactor{i}'])], axis=1)
complex_b_genes = interactions.apply(lambda row: [row[f'interactor{i}'] for i in range(row['num_interactors_a'] + 1, row['num_interactors_a'] + row['num_interactors_b'] + 1) if pd.notna(row[f'interactor{i}'])], axis=1)
interactions_all_genes = complex_a_genes + complex_b_genes

# remove duplicates from each element of interactions_all_genes
interactions_all_genes = interactions_all_genes.apply(lambda x: list(set(x)))

# --------- interactions ---------

genes = WGCNA.datExpr.var

# Only keep interactions where all genes are in the WGCNA data
print(f'Number of interactions before filtering: {len(interactions)}')
interactions = interactions[interactions_all_genes.apply(lambda x: all([gene in genes.index for gene in x]))] 
print("Keeping interactions with all genes present: ", len(interactions))

# Remove interactions with only one gene
interactions = interactions[interactions_all_genes.apply(lambda x: len(x) > 1)]
print("Keeping interactions with more than one gene: ", len(interactions))

result = interactions.copy()

result["same_module"] = False
result["module"] = None
result["mean_adj"] = None

for i, row in interactions.iloc[:, 3:].iterrows():
    all_genes = [gene for gene in row if pd.notna(gene)]

    if WGCNA.datExpr.var.loc[all_genes, "moduleLabels"].nunique() == 1:
        result.loc[i, "same_module"] = True
        result.loc[i, "module"] = genes.loc[all_genes[0], "moduleLabels"]
    else:
        result.loc[i, "same_module"] = False
        result.loc[i, "module"] = None
    # Add adjacency to interactions
    adj = WGCNA.adjacency.loc[all_genes, all_genes]
    idxs_x, idxs_y = np.triu_indices(adj.shape[0], 1)
    adj = adj.values[idxs_x, idxs_y].mean()
    result.loc[i, "mean_adj"] = adj
    
print("Showing the first 10 interactions")
print(result.head(10))

# Save the interactions in the same path as the WGCNA object
print("Saving interactions to {}".format(os.path.join(output_path, "interactions.csv")))
result.to_csv(os.path.join(output_path, "interactions.csv"), index=False)

print("Done: interactions")

# ----- Occurrences ------

# get all genes that appear in at least one interaction
all_interactors = set()
for interaction in interactions_all_genes:
    all_interactors.update(interaction)

all_interactors = all_interactors.intersection(genes.index)
all_interactors = list(all_interactors)

# Count occurrences of interactors in modules
interactors_counts = pd.DataFrame(0, index=all_interactors, columns=genes["moduleLabels"].unique())

for gene in all_interactors:
    module = genes.loc[gene, "moduleLabels"]
    interactors_counts.loc[gene, module] += 1

print("Saving interactors occurrences to: ", os.path.join(output_path, "interactors_occurrences.csv"))
interactors_counts.to_csv(os.path.join(output_path, "interactors_occurrences.csv"))

# ----- Module info ------

# rank modules by size and save number of LR_interactions in each module
print("Ranking modules by size and saving number of LR interactions in each module")
# number of modules
n_modules = len(WGCNA.datExpr.var["moduleLabels"].unique())
# make dataframe with module sizes
module_sizes = WGCNA.datExpr.var["moduleLabels"].value_counts().to_dict()
module_df = pd.DataFrame(data=module_sizes.items(), columns=["module", "size"])
# sort by size, bigger first
module_df = module_df.sort_values(by="size", ascending=False)
module_df["rank"] = range(1, module_df.shape[0] + 1)
# add number of LR interactions in each module
module_df["n_interactions"] = None

for i, row in module_df.iterrows():
    n_interactions = result[result["module"] == row["module"]]["same_module"].sum()
    module_df.loc[i, "n_interactions"] = n_interactions
    
module_df.to_csv(os.path.join(output_path, "module_info.csv"), index=False)

print("Done: Module info")

# ----- General info ------
# If interactions_info.txt exists, overwrite it
print("interactions_info.txt already exists, overwriting")
if os.path.exists(os.path.join(output_path, "interactions_info.txt")):
    os.remove(os.path.join(output_path, "interactions_info.txt"))

print("Writing general info")

# append new lines, don't overwrite
with open(os.path.join(output_path, "interactions_info.txt"), "w") as f:
    f.write("LR_resource: {}\n".format(args.lr_resource))
    f.write("total_interactions: {}\n".format(result.shape[0]))
    f.write("n_interactions_same_module: {}\n".format(result["same_module"].sum()))

# ----- Done ------

print("Done: interactions.py")
