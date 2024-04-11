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

print("Variance stabilizing transformation")
# Perform variance stabilizing transformation
dds.vst(use_design=False)

# Remove dummy design factors
if design_factors == "dummy_deseq2_condition": 
    adata.obs.drop(columns=design_factors, inplace=True)

# Get results back to adata object
adata.X = dds.layers["vst_counts"].copy()
adata.layers["deseq2_vst_counts"] = dds.layers["vst_counts"].copy()
adata.layers["deseq2_norm_counts"] = dds.layers["normed_counts"].copy()

# Overwrite adata object
print("Overwriting adata object")
adata.write(args.input, compression="gzip")

print("Done: deseq_norm.py")

