import numpy as np
import sys
import scanpy as sc
import pandas as pd
import pingouin as pg
import os
import argparse
import ast

parser = argparse.ArgumentParser(description='Run WGCNA pipeline')

parser.add_argument('-i', '--input', type=str, help='path to input adata file')

args = parser.parse_args()

filename = args.input

# name is the last part of the filename after the last /
name = filename.split("/")[-1].split(".")[0]

# output_path is the path to the folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("name: ", name)
print("output_path: ", output_path)

## --------- CORRELATIONS ------------
#
#print("Reading adata file: {}".format(filename))
#adata = sc.read_h5ad(filename)
#
## make expression dataframe
#geneExp = pd.DataFrame(
#    data=adata.X,
#    index=adata.obs_names,
#    columns=adata.var_names
#)
#
## Pearson correlation
#corr_matrix = geneExp.corr()
#outfile = os.path.join(output_path, 'correlation.csv.gz')
#corr_matrix.to_csv(outfile, compression='gzip')
#
## Partial correlation using :contentReference[oaicite:0]{index=0}
#partial_corr_matrix = pg.pcorr(geneExp)
#outfile = os.path.join(output_path, 'partial_correlation.csv.gz')
#partial_corr_matrix.to_csv(outfile, compression='gzip')
#
#print("Correlation matrices saved")
#
# --------- ADD PCORR TO INTERACTIONS ------------

# corr and pcorr are already computed, just load pcorr
partial_corr_file = os.path.join(output_path, 'partial_correlation.csv.gz')
partial_corr_matrix = pd.read_csv(partial_corr_file, index_col=0)

print("\nAdding pcorr to selected interaction results")

interactions_resources = {
    'ccc_lr_pairs': '/home/lnemati/pathway_crosstalk/data/interactions/ccc_lr_pairs.csv',
    'intact_direct': '/home/lnemati/pathway_crosstalk/data/interactions/intact_direct.csv',
    'intact_physical': '/home/lnemati/pathway_crosstalk/data/interactions/intact_physical.csv',
    'intact_association': '/home/lnemati/pathway_crosstalk/data/interactions/intact_association.csv',
}

for name, _ in interactions_resources.items():

    existing_file = os.path.join(output_path, f"interactions/{name}.csv")

    if not os.path.exists(existing_file):
        print(f"Skipping {name}: scored file not found")
        continue

    print(f"Processing {name}")

    df = pd.read_csv(existing_file)

    if "all_genes" not in df.columns:
        print("  skipped (no all_genes column)")
        continue

    df["all_genes"] = df["all_genes"].apply(
        lambda x: ast.literal_eval(x) if isinstance(x, str) else x
    )

    pcorr_values = []

    for genes in df["all_genes"]:

        if not isinstance(genes, (list, tuple)) or len(genes) < 2:
            pcorr_values.append(np.nan)
            continue

        g1, g2 = genes[0], genes[1]

        if g1 in partial_corr_matrix.index and g2 in partial_corr_matrix.columns:
            pcorr_values.append(partial_corr_matrix.loc[g1, g2])
        else:
            pcorr_values.append(np.nan)

    # add pcorr column WITHOUT modifying other data
    df["pcorr"] = pcorr_values

    new_file = os.path.join(output_path, f"interactions/{name}_with_pcorr.csv")
    df.to_csv(new_file, index=False)

    print(f"  saved → {new_file}")

print("Done: partial_corr.py")
