import numpy as np
import sys
import scanpy as sc
import pandas as pd
import os
import argparse
from pydeseq2.dds import DeseqDataSet, DefaultInference

# read filename from command line
#filename = os.environ.get("WGCNA_FILE_PATH")

parser = argparse.ArgumentParser(description='Normalize adata using DESeq2')

parser.add_argument('-i', '--input', type=str, help='path to input adata file')
parser.add_argument('-n', '--ncpus', type=int, help='number of cpus to use')

args = parser.parse_args()

filename = args.input
adata = sc.read_h5ad(filename)
n_cpus = args.ncpus

# --------- DeSeq2 -----------

print("Starting DeSeq2 normalization")

adata.X = adata.layers["raw_counts"].todense()

design_factors = "condition" if len(adata.obs["condition"].unique()) > 1 else None

if design_factors is None:
    print("Cannot run Deseq2 on a single condition, creating dummy design factors")
    # Create dummy design factors that alternate between A and B
    adata.obs["dummy_deseq2_condition"] = ["A" if i % 2 == 0 else "B" for i in range(len(adata.obs))]
    design_factors = "dummy_deseq2_condition"

print("Initializing Deseq2")
# Initialize DeseqDataSet
inference = DefaultInference(n_cpus=n_cpus)

dds = DeseqDataSet(
    counts=pd.DataFrame(data=adata.X.astype(int), index=adata.obs_names, columns=adata.var_names),
    metadata=adata.obs,
    design_factors=design_factors,
    refit_cooks=True,
    inference=inference
)

print("Running Deseq2 normalization")
# Run DESeq2 normalization
dds.fit_size_factors()

# Get results back to adata object
adata.layers["deseq2_norm_counts"] = dds.layers["normed_counts"].copy()
adata.X = adata.layers["deseq2_norm_counts"].copy()

print("Variance stabilizing transformation")
# Perform variance stabilizing transformation
dds.vst(use_design=False)
adata.layers["deseq2_vst_counts"] = dds.layers["vst_counts"].copy()

# Remove dummy design factors
if design_factors == "dummy_deseq2_condition": 
    adata.obs.drop(columns=design_factors, inplace=True)

# --------- PCA Outliers -----------

# Calculate PCA
sc.pp.pca(adata, n_comps=10, use_highly_variable=False) 

# Calculate centroid of the PCA using median
centroid = np.median(adata.obsm["X_pca"], axis=0)

# Calculate distance from centroid
adata.obs["pca_distance"] = np.linalg.norm(adata.obsm["X_pca"] - centroid, axis=1)

# Remove samples that are more than 5 stds away from the median
std = np.std(adata.obs["pca_distance"])

adata.obs['pca_outlier'] = adata.obs["pca_distance"] > 5 * std
print("Found {} outliers".format(adata.obs['pca_outlier'].sum()))

# If pca_outliers are more than 5% of the samples, print a warning to stderr
if adata.obs['pca_outlier'].sum() > 0.05 * len(adata.obs):
    print("WARNING: More than 5% of the samples are outliers", file=sys.stderr)

# Remove outliers
adata = adata[~adata.obs['pca_outlier']]
adata.obs.drop(columns=["pca_distance", "pca_outlier"], inplace=True)

# Remove pca from obsm
adata.obsm.pop("X_pca")
# Overwrite adata object

# ----- Save adata object ------

# Ser adata.X to vst
adata.X = adata.layers["deseq2_vst_counts"].copy()

print("Overwriting adata object")
adata.write(args.input, compression="gzip")

print("Done: deseq_norm.py")

