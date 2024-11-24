import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Get values for gpcrs, gproteins and ligands')

parser.add_argument('-i', '--input', type=str, help='path to input wgcna file')

args = parser.parse_args()

filename = args.input
WGCNA = PyWGCNA.readWGCNA(filename)

output_path = "/".join(filename.split("/")[:-1]) + "/"
condition = filename.split("/")[-3]

print("filename", filename)
print("condition", condition)
print("output_path", output_path)

# Read GPCRs
print("Reading GPCRs")
gpcrs = pd.read_csv('/home/fraimondi/BIOINFO1_BK/francesco/GPCR/gpcrs_genes.txt')['Symbol'].values

# Read ligands
print("Reading ligands")
lr = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc_lr_pairs.csv')
lr = lr[lr.receptor.isin(gpcrs)]

# Read G-Proteins coupling
print("Reading G-Proteins")
gprots_couplings = pd.read_csv('/home/lnemati/resources/gpcr_couplings/gpcr_couplings.csv')
gprots_couplings = gprots_couplings.query('type != "precogx_only"')
gprots_couplings.set_index('gpcr', inplace=True, drop=False)
gprots_couplings.index.name = None
gprots_couplings = gprots_couplings.groupby('gpcr')['gprot'].apply(list)
gprots = gprots_couplings.values.sum()

# Group all together
all_gpcrs_genes = list(set(gpcrs) | set(gprots) | set(lr.ligand.values) | set(lr.receptor.values))
all_gpcrs_genes = WGCNA.datExpr.var.index.intersection(all_gpcrs_genes)
print("Number of genes", len(all_gpcrs_genes))

# Get correlation
corr = pd.DataFrame(np.corrcoef(WGCNA.datExpr.X.T), index=WGCNA.datExpr.var.index, columns=WGCNA.datExpr.var.index) 
corr = corr.loc[all_gpcrs_genes, all_gpcrs_genes]

# Normalize adjacency so that maximum degree is 1
adj = WGCNA.adjacency.loc[all_gpcrs_genes, all_gpcrs_genes]
adj = adj / adj.sum().max()

# Make a dataframe with all  pairs
df = pd.DataFrame(columns=['gpcr', 'interactor', 'condition', 'interactor_type', 'value_type', 'value'])

def create_df(matrix, value_type):
    df = pd.DataFrame(columns=['gpcr', 'interactor', 'condition', 'interactor_type', 'value_type', 'value'])

    # Read and add GPCR-G-Protein
    valid_gpcrs = matrix.index.intersection(gpcrs).intersection(gprots_couplings.index)
    for gpcr in valid_gpcrs:
        couplings = matrix.index.intersection(gprots_couplings.loc[gpcr])
        values = matrix.loc[gpcr, couplings].to_frame(name = 'value')
        values['gpcr'] = gpcr
        values['interactor'] = values.index
        values['condition'] = condition
        values['interactor_type'] = 'G-Protein'
        values['value_type'] = value_type
        df = pd.concat([df, values])

    # Read and add GPCR-Ligand 
    valid_gpcrs = matrix.index.intersection(gpcrs).intersection(lr.receptor.unique())
    for gpcr in valid_gpcrs:
        ligands = matrix.index.intersection(lr[lr['receptor'] == gpcr].ligand.unique())
        values = matrix.loc[gpcr, ligands].to_frame(name = 'value')
        values['gpcr'] = gpcr
        values['interactor'] = values.index
        values['condition'] = condition
        values['interactor_type'] = 'Ligand'
        values['value_type'] = value_type
        df = pd.concat([df, values])

    df = df.reset_index(drop=True)

    return df

df = pd.DataFrame(columns=['gpcr', 'interactor', 'condition', 'interactor_type', 'value_type', 'value'])
# Create the dataframe
df = pd.concat([df, create_df(corr, 'correlation'), create_df(adj, 'adjacency')])

# Save the dataframe
df.to_csv(os.path.join(output_path, 'gpcrs.csv'), index=False)

print("Done: gpcrs.py")
