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

# Get adjacency matrix and normalize it so that the maximum degree is 1
adj = wgcna.adjacency
degree = adj.sum(axis=1)
adj = adj / degree.max()

# Get correlation matrix
corr = np.corrcoef(wgcna.datExpr.X.T)
corr = pd.DataFrame(corr, index=wgcna_genes.index, columns=wgcna_genes.index)

# output_path is the path to the folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("output_path", output_path)
print("filename", filename)

interactions_resources = {
    #'ccc':'/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv',
    #'pairs_of_interactions':'/home/lnemati/pathway_crosstalk/data/interactions/pairs_of_interactions.csv',
    #'all_ccc_gene_pairs':'/home/lnemati/pathway_crosstalk/data/interactions/all_ccc_gene_pairs.csv',
    'ccc_lr_pairs': '/home/lnemati/pathway_crosstalk/data/interactions/ccc_lr_pairs.csv',
    'intact_direct':'/home/lnemati/pathway_crosstalk/data/interactions/intact_direct.csv',
    'intact_physical':'/home/lnemati/pathway_crosstalk/data/interactions/intact_physical.csv',
    'intact_association':'/home/lnemati/pathway_crosstalk/data/interactions/intact_association.csv',
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

    # TODO: define here the different scores
    result["module"] = None
    result["same_module"] = 0
    result['corr'] = None
    result['adj'] = None
    #result["min_adj"] = None
    #result["mean_adj"] = None
    result['kme_corr'] = None
    #result["min_kme_corr"] = None
    #result["mean_kme_corr"] = None
    result["missing_genes"] = False

    for i, row in interactions.iterrows():
        all_genes = row["all_genes"]
        # Check if all genes are in the WGCNA object
        genes_in_wgcna = [gene for gene in all_genes if gene in wgcna_genes.index]
        
        if genes_in_wgcna != all_genes:
            # Check for missing genes
            result.loc[i, "missing_genes"] = True
        else: 
            # If all genes are present, calculate the scores
            result.loc[i, 'corr'] = corr.loc[all_genes[0], all_genes[1]]
            result.loc[i, 'adj'] = adj.loc[all_genes[0], all_genes[1]]
            kme1 = wgcna.signedKME.loc[all_genes[0]]
            kme2 = wgcna.signedKME.loc[all_genes[1]]
            result.loc[i, 'kme_corr'] = kme1.corr(kme2)
            if wgcna.datExpr.var.loc[genes_in_wgcna, "moduleLabels"].nunique() == 1:
                # If all genes are present and in the same module add module
                result.loc[i, "same_module"] = 1
                result.loc[i, "module"] = wgcna_genes.loc[all_genes[0], "moduleLabels"]
        
        # This made sense when considering interactions, now that we use pairs, we don't need it

        # Add adjacency to interactions, use 0 value to fill for genes not in the WGCNA object
        #adj = wgcna.adjacency.reindex(index=all_genes, columns=all_genes, fill_value=0)
        #idxs_x, idxs_y = np.triu_indices(adj.shape[0], 1)
        #result.loc[i, "mean_adj"] = np.mean(adj.values[idxs_x, idxs_y])
        #result.loc[i, "min_adj"] = np.min(adj.values[idxs_x, idxs_y])
        
        # Add KME correlation to interactions, use 0 value to fill for genes not in the WGCNA object
        #kme = wgcna.signedKME.loc[genes_in_wgcna]
        #corr = kme.T.corr().reindex(index=all_genes, columns=all_genes, fill_value=0)
        #idxs_x, idxs_y = np.triu_indices(corr.shape[0], 1)
        #result.loc[i, "mean_kme_corr"] = np.mean(corr.values[idxs_x, idxs_y])
        #result.loc[i, "min_kme_corr"] = np.min(corr.values[idxs_x, idxs_y])
    
    # Show results
    print("Showing the first 10 interactions")
    print(result.head(10))

    # Create 'interactions' folder if it doesn't exist
    os.makedirs(os.path.join(output_path, "interactions"), exist_ok=True)

    # Save results
    result.to_csv(os.path.join(output_path, f"interactions/{name}.csv"), index=False)

    print("Done: score_interactions.py")
