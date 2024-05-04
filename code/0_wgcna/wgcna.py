import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
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

WGCNA = PyWGCNA.WGCNA(
    name="wgcna_" + name, 
    species='homo sapiens', 
    #anndata=adata, 
    geneExp=geneExp,
    TPMcutoff=0,
    RsquaredCut=0.85,
    networkType="signed hybrid",
    minModuleSize=50,
    powers=list(range(1,11)) + list(range(12,31))[::2],
    outputPath=output_path,
    save=True
)

# add sample and gene info
WGCNA.updateSampleInfo(adata.obs)
WGCNA.updateGeneInfo(adata.var)

WGCNA.preprocess()
WGCNA.findModules()
WGCNA.saveWGCNA()

print("Done: wgcna.py")
