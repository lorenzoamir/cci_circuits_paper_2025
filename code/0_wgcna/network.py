import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import powerlaw
import os
import igraph as ig
from scipy.sparse import csr_matrix
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

sft_power = WGCNA.power
adj = WGCNA.adjacency # Pandas DataFrame

# Remove self-loops by filling diagonal with 0
#np.fill_diagonal(adj.values, 0)

degree_df = pd.DataFrame(index=adj.index, columns=["degree", "intramodular_degree", "module"])
degree_df["degree"] = adj.sum(axis=0)

# For each module, get the sum of adjacency matrix for all genes in the module
# Subset both rows and columns so that you only get the adjacency matrix for the module
for module, df in adj.groupby(WGCNA.datExpr.var['moduleLabels'], axis=0):
    degree_df.loc[df.index, "module"] = module
    degree_df.loc[df.index, "intramodular_degree"] = df.loc[df.index, df.index].sum(axis=0)

# Fit degree distribution to powerlaw
results = powerlaw.Fit(degree_df["degree"].values)
powerlaw_alpha = results.power_law.alpha
powerlaw_xmin = results.power_law.xmin
powerlaw_xmax = results.power_law.xmax
powerlaw_sigma = results.power_law.sigma

# Create an igraph graph from the adjacency matrix
# Remove weak interactions that sum up to 10% of the total weight

# Flatten the matrix
flat_adj = adj.values
flat_adj = flat_adj[flat_adj > 0]
flat_adj = np.sort(flat_adj)

# Get total and cutoff values
total_weight = np.sum(flat_adj)
cutoff = 0.10 * total_weight

# Get weakest link to keep
cum_sum = np.cumsum(flat_adj)
num_to_remove = np.searchsorted(cum_sum, cutoff, side='right')
weakest_link = flat_adj[num_to_remove]

g = ig.Graph.Weighted_Adjacency(csr_matrix(np.where(adj < weakest_link, 0, adj)), mode=ig.ADJ_UNDIRECTED)

# Distance = 1 - adj
d = 1 - np.array(g.es['weight'])

# Distance must be strictly positive, so set 0 to min/10
d[d == 0] = np.min(d[d > 0]) / 10

# Set distance attribute
g.es['distance'] = d

# Calculate network attributes
clustering_global = g.transitivity_undirected()
modularity = g.modularity(degree_df["module"], weights=g.es['weight'])
avg_path_length = g.average_path_length(weights=g.es['distance'])
diameter = g.diameter(directed=False, weights=g.es['distance'])

# Calculate node attributes
degree_df["pagerank"] = g.pagerank(weights=g.es['weight'])
degree_df["clustering_local"] = g.transitivity_local_undirected(weights=g.es['weight'])
degree_df["closeness"] = g.closeness(weights=g.es['distance'])
degree_df["betweenness"] = g.betweenness(weights=g.es['distance'])

# save degree_df as csv
degree_df.to_csv(os.path.join(output_path, "degree_df.csv.gz"))

# Maybe add assortativity?


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
    
    f.write("sft_power: {}\n".format(sft_power))
    
    f.write("n_pairs_all_genes: {}\n".format(total_pairs))
    f.write("n_pairs_same_module_all_genes: {}\n".format(n_same_module))
    
    f.write("avg_degree: {}\n".format(degree_df["degree"].mean()))
    f.write("avg_intramodular_degree: {}\n".format(degree_df["intramodular_degree"].mean()))
    f.write("var_degree: {}\n".format(degree_df["degree"].var()))
    f.write("var_intramodular_degree: {}\n".format(degree_df["intramodular_degree"].var()))

    f.write("powerlaw_alpha: {}\n".format(powerlaw_alpha))
    f.write("powerlaw_xmin: {}\n".format(powerlaw_xmin))
    f.write("powerlaw_xmax: {}\n".format(powerlaw_xmax))
    f.write("powerlaw_sigma: {}\n".format(powerlaw_sigma))

    f.write("avg_path_length: {}\n".format(avg_path_length))
    f.write("diameter: {}\n".format(diameter))
    f.write("clustering_global: {}\n".format(clustering_global))
    f.write("modularity: {}\n".format(modularity))

# ----- Done ------

print("Done: network.py")
