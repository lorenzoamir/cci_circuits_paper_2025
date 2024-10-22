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
wgcna = PyWGCNA.readWGCNA(filename)
wgcna.CalculateSignedKME()
wgcna_genes = wgcna.datExpr.var

# output_path is the path to the folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("output_path", output_path)
print("filename", filename)

interactions_resources = {
    'ccc':'/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv',
    'intact_direct':'/home/lnemati/pathway_crosstalk/data/interactions/intact_direct.csv',
}

for name, interaction_path in interactions_resources.items():
    print(f"Loading {name} interactions from {interaction_path}")
    interactions = pd.read_csv(interaction_path)
    print(f"Number of interactions: {len(interactions)}")

    # Convert all_genes to list
    interactions['all_genes'] = interactions['all_genes'].apply(ast.literal_eval)
    print('First 3 interactions:')
    print(interactions.head(3))

    result = interactions.copy()

    result["module"] = None
    result["same_module"] = None
    result["min_adj"] = None
    result["median_adj"] = None
    result["min_kme_corr"] = None
    result["median_kme_corr"] = None

    for i, row in interactions.iterrows():
        all_genes = row["all_genes"]
        # Check if all genes are in the WGCNA object
        genes_in_wgcna = [gene for gene in all_genes if gene in wgcna_genes.index]

        if wgcna.datExpr.var.loc[genes_in_wgcna, "moduleLabels"].nunique() == 1 and genes_in_wgcna == all_genes:
            result.loc[i, "same_module"] = 1
            result.loc[i, "module"] = wgcna_genes.loc[all_genes[0], "moduleLabels"]
        else:
            result.loc[i, "same_module"] = 0
            result.loc[i, "module"] = None

        # Add adjacency to interactions, use 0 value to fill for genes not in the WGCNA object
        adj = wgcna.adjacency.reindex(index=all_genes, columns=all_genes, fill_value=0)
        idxs_x, idxs_y = np.triu_indices(adj.shape[0], 1)
        result.loc[i, "median_adj"] = np.median(adj.values[idxs_x, idxs_y])
        result.loc[i, "min_adj"] = np.min(adj.values[idxs_x, idxs_y])

        # Add KME correlation to interactions, use 0 value to fill for genes not in the WGCNA object
        kme = wgcna.signedKME.loc[genes_in_wgcna]
        corr = kme.T.corr().reindex(index=all_genes, columns=all_genes, fill_value=0)
        idxs_x, idxs_y = np.triu_indices(corr.shape[0], 1)
        result.loc[i, "median_kme_corr"] = np.median(corr.values[idxs_x, idxs_y])
        result.loc[i, "min_kme_corr"] = np.min(corr.values[idxs_x, idxs_y])
    
    # Show results
    print("Showing the first 10 interactions")
    print(result.head(10))

    # Create 'interactions' folder if it doesn't exist
    os.makedirs(os.path.join(output_path, "interactions"), exist_ok=True)

    # Save results
    result.to_csv(os.path.join(output_path, f"interactions/{name}.csv"), index=False)

    print("Done: score_interactions.py")
