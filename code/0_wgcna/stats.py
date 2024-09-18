import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import PyWGCNA
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import mannwhitneyu
from itertools import combinations, product
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Perform statistical analysis on WGCNA results')

parser.add_argument('-i', '--inputfile', dest='inputfile', help='Input WGCNA object file', required=True)

args = parser.parse_args()

print("inputfile: " + args.inputfile)

# Create output directory
roc_dir = '/home/lnemati/pathway_crosstalk/results/roc'
os.makedirs(roc_dir, exist_ok=True)
os.makedirs(os.path.join(roc_dir, "tumor"), exist_ok=True)
os.makedirs(os.path.join(roc_dir, "normal"), exist_ok=True)

# read interactions
workdir = os.path.dirname(args.inputfile)
figdir = os.path.join(workdir, "figures")

# file to write results to
stats_file = os.path.join(workdir, "stats.txt")

# create it if it does not exist, overwrite it if it does
with open(stats_file, "w") as f:
    f.write("")

# Function that generates and saves the ROC curve
def generate_roc_curve(
    data,
    target_col,
    feature_col,
    file,
    name=None,
):
    # ROC curve using all gene pairs
    fpr, tpr, ths = roc_curve(data[target_col], data[feature_col])
    auroc = roc_auc_score(data[target_col], data[feature_col])
    
    with open(file, "a") as f:
        f.write(f"auroc: " + str(auroc) + "\n")

    return fpr, tpr, ths

def update_diff_complex(
    all_pairs,
    genes_series_a,
    genes_series_b
):
    for genes_a, genes_b in zip(genes_series_a, genes_series_b):
        for pair in product(genes_a, genes_b):
            gene_a, gene_b = pair
            if pair in all_pairs.index:
                all_pairs.loc[pair, "interaction"] = 1
            elif (gene_b, gene_a) in all_pairs.index:
                all_pairs.loc[(gene_b, gene_a), "interaction"] = 1
    return

print("Reading WGCNA ojbects")
WGCNA = PyWGCNA.readWGCNA(args.inputfile)
print("Reading interactions.csv")
# read interactions.csv 
interactions = pd.read_csv(os.path.join(workdir, "interactions.csv"))
# make dict with all genes in each module
modules_t = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}

# Extract genes for complex A and B for each row
complex_a_genes = interactions.apply(lambda row: list(set([row[f'interactor{i}'] for i in range(1, row['num_interactors_a'] + 1) if pd.notna(row[f'interactor{i}'])])), axis=1)
complex_b_genes = interactions.apply(lambda row: list(set([row[f'interactor{i}'] for i in range(row['num_interactors_a'] + 1, row['num_interactors_a'] + row['num_interactors_b'] + 1) if pd.notna(row[f'interactor{i}'])])), axis=1)
interactions_all_genes = complex_a_genes + complex_b_genes

# make all_interacting_genes a list containing all unique genes in the interactions
all_interacting_genes = list(set().union(*interactions_all_genes))
print("Number of interacting genes: ", len(all_interacting_genes))

# ----- ROC all gene pairs -----

print("Generating ROC curves")

adj = WGCNA.adjacency
print("Adjacency has shape: ", adj.shape)

# flatten matrix, remove diagonal and duplicated values
print("Flattening matrix")
all_pairs = pd.DataFrame(
    adj.where(
        np.tri(
            adj.shape[0],
            dtype=bool,
            k=-1
        ),
        np.nan
    ).stack(dropna=True), columns=["adj"]
)
print("all_pairs has now shape: ", all_pairs.shape)
print(all_pairs.head())

all_pairs['interaction'] = 0

update_diff_complex(all_pairs, complex_a_genes, complex_b_genes)

print("all_pairs has now shape: ", all_pairs.shape)
print("Total number of interacting pairs: ", all_pairs["interaction"].sum())

print(all_pairs.head())

# find nans in all_pairs
print("Number of nans in all_pairs: ", all_pairs.isna().sum())

name = WGCNA.name
if "/tumor/" in args.inputfile:
    condition = "tumor"
if "/normal/" in args.inputfile:
    condition = "normal"

# Diff complex interactions
fpr, tpr, ths = generate_roc_curve(
    data=all_pairs,
    target_col="interaction",
    feature_col="adj",
    file=stats_file,
    name=name
)

# save fpr, tpr, ths
print("Saving ROC curve")
output_dir = os.path.join(roc_dir, condition, name)
os.makedirs(output_dir, exist_ok=True)
print("Saving to: " + output_dir)
np.save(os.path.join(output_dir, "fpr.npy"), fpr)
np.save(os.path.join(output_dir, "tpr.npy"), tpr)
np.save(os.path.join(output_dir, "ths.npy"), ths)
print()

# ------ Test VS All -------
print("Performing test")

# Perform test: do interacting pairs have higher connectivity than other pairs?
U_all, p_all = mannwhitneyu(
    all_pairs.loc[all_pairs['interaction'] == True, "adj"],
    all_pairs.loc[all_pairs['interaction'] == False, "adj"],
    alternative="greater"
)
print("Writing result to: " + stats_file)

# create it if it does not exist, overwrite it if it does
with open(stats_file, "a") as f:
    f.write("mannwhitneyu_U_all: " + str(U_all) + "\n")
    f.write("mannwhitneyu_p_all: " + str(p_all) + "\n")
    
print("Done: stats.py") 
