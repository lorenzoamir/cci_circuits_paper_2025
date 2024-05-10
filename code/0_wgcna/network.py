import numpy as np
import sys
import scanpy as sc
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

# --------- Check corresponding normal tissue -----------
print("Checking corresponding normal tissue")
normal_tissue_files = None

if "/tumor/" in filename:
    print("Checking corresponding normal tissue")
    # Get the corresponding normal tissue, strip everything after /tumor/
    normal_tissue = filename.split("/tumor/")[0] + "/normal/"
    # Search for .h5ad files in the normal tissue folder
    # and save their absolute path
    normal_tissue_files = [os.path.join(normal_tissue, f) for f in os.listdir(normal_tissue) if f.endswith(".h5ad")]
    # replace .h5ad with .p
    normal_tissue_files = [f.replace(".h5ad", ".p") for f in normal_tissue_files]
    # Check if there is only one corresponding normal tissue
    if len(normal_tissue_files) == 1:
        print("Found corresponding normal tissue: {}".format(normal_tissue_files[0]))
        # Writing normal tissue location in same folder as tumor tissue
        with open(os.path.join(output_path, "corresponding_normal_file.txt"), "w") as f:
            f.write(normal_tissue_files[0])
    elif len(normal_tissue_files) == 0:
        print("No corresponding normal tissue found")
        normal_tissue_files = None
    else:
        print("Multiple corresponding normal tissues found: {}".format(normal_tissue_files))
        normal_tissue_files = None

print("Done: Checking corresponding normal tissue")

# ------ TOM & adjacency ------

print("Saving TOM and adjacency")

# save TOM matrix as csv
WGCNA.TOM.to_csv(os.path.join(output_path, "tom.csv.gz"))
# save adjacency matrix as csv
WGCNA.adjacency.to_csv(os.path.join(output_path, "adjacency.csv.gz"))

print("Done: TOM and adjacency")

# -------- Same Module ---------

print("Checking which gene pairs are in the same module")

same_module = WGCNA.datExpr.var['moduleLabels'].values[:, None] == WGCNA.datExpr.var['moduleLabels'].values
same_module = pd.DataFrame(data=same_module, index=WGCNA.datExpr.var.index, columns=WGCNA.datExpr.var.index)

# save same_module as csv
same_module.to_csv(os.path.join(output_path, "same_module.csv.gz"))

# Get idxs corresponding to upper triangle (diagonal excluded)
idxs_x, idxs_y = np.triu_indices(same_module.shape[0], 1)

same_module_vals = same_module.values[idxs_x, idxs_y]

total_pairs = same_module_vals.shape[0]
n_same_module = same_module_vals.sum()

# ----- Network ------

adj = WGCNA.adjacency

degree_df = pd.DataFrame(index=adj.index, columns=["degree", "intramodular_degree", "module"])
degree_df["degree"] = adj.sum(axis=0)

# For each module, get the sum of adjacency matrix for all genes in the module
# Subset both rows and columns so that you only get the adjacency matrix for the module
for module, df in adj.groupby(WGCNA.datExpr.var['moduleLabels'], axis=0):
    degree_df.loc[df.index, "module"] = module
    degree_df.loc[df.index, "intramodular_degree"] = df.loc[df.index, df.index].sum(axis=0)

# save degree_df as csv
degree_df.to_csv(os.path.join(output_path, "degree_df.csv.gz"))

# ----- General info ------

# If general_info.txt exists, overwrite it
print("general_info.txt already exists, overwriting")
if os.path.exists(os.path.join(output_path, "general_info.txt")):
    os.remove(os.path.join(output_path, "general_info.txt"))

print("Writing general info")
with open(os.path.join(output_path, "general_info.txt"), "w") as f:
    f.write("name: {}\n".format("wgcna_" + category))
    f.write("wgcna_file: {}\n".format(os.path.join(output_path, "wgcna_" + category + ".p")))
    if normal_tissue_files:
        if len(normal_tissue_files) == 1:
            f.write("corresponding_normal_tissue: {}\n".format(normal_tissue_files[0]))

    f.write("n_modules: {}\n".format(WGCNA.datExpr.var['moduleLabels'].nunique()))
    f.write("module_sizes: {}\n".format(WGCNA.datExpr.var['moduleLabels'].value_counts().to_dict()))

    f.write("n_genes: {}\n".format(WGCNA.datExpr.shape[1]))
    f.write("n_samples: {}\n".format(WGCNA.datExpr.shape[0]))

    f.write("n_pairs_all_genes: {}\n".format(total_pairs))
    f.write("n_pairs_same_module_all_genes: {}\n".format(n_same_module))

# ----- Done ------

print("Done: network.py")
