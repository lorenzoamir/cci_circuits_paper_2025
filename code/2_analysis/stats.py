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
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumor', type=str, help='Input tumor WGCNA')
parser.add_argument('--normal', type=str, help='Input normal WGCNA')

args = parser.parse_args()

print("normal file:", args.normal)
print("tumor file:", args.tumor)

# read general_info.txt in the tumor directory
info_dir = args.tumor.split("/")[0:-1]
info_file = os.path.join("/".join(info_dir), "general_info.txt")
print("Trying to read: " + info_file)

# read interactions

tumor_dir = os.path.dirname(args.tumor)
normal_dir = os.path.dirname(args.normal)

tumor_fig_dir = os.path.join(os.path.dirname(args.tumor), "figures")
normal_fig_dir = os.path.join(os.path.dirname(args.normal), "figures")

comparison_file = os.path.join(tumor_dir, "comparison.txt")

# Create comparison.txt file in the tumor directory if it does not exist
with open(comparison_file, "w") as f:
    f.write("Tumor directory: " + tumor_dir + "\n")
    f.write("Normal directory: " + normal_dir + "\n")

print("Tumor directory: " + tumor_dir)
print("Normal directory: " + normal_dir)
print("Tumor figures will be saved to: " + tumor_fig_dir)
print("Normal figures will be saved to: " + normal_fig_dir)
print("Comparison file: " + comparison_file)

# Plot heatmap of TOM for tumor and normal comparing all genes with ligand and receptors only
fs = 18
cmap='viridis'

# Create figures for the ROC curves
roc_fig, roc_ax = plt.subplots(1, 1, figsize=(10,10))

# Function that generates and saves the ROC curve
def generate_roc_curve(
    data,
    target_col,
    feature_col,
    condition,
    comparison_file,
    ax,
    name=None,
):
    # ROC curve using all gene pairs
    fpr, tpr, _ = roc_curve(data[target_col], data[feature_col])
    auroc = roc_auc_score(data[target_col], data[feature_col])
    
    if condition.lower() == "tumor":
        with open(comparison_file, "a") as f:
            f.write(f"auroc_tumor_{target_col}: " + str(auroc) + "\n")
    elif condition.lower() == "normal":
        with open(comparison_file, "a") as f:
            f.write(f"auroc_normal_{target_col}: " + str(auroc) + "\n")

    if condition == "tumor":
        color = 'C1' 
        label = f'TCGA {target_col} (AUC = {auroc:.2f})'
    elif condition == "normal":
        color = 'C0'
        label = f'GTEx {target_col} (AUC = {auroc:.2f})'

    # add title to the plot if name is provided
    if name:
        ax.set_title(name + " ROC curve", fontsize=fs)
    ax.plot([0, 1], [0, 1], color='k', lw=2)
    ax.plot(fpr, tpr, color=color, lw=2, label=label, )
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend()
    
    return

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


for condition in ["normal", "tumor"]:

    if condition == "tumor":
        print("===== Tumor =====\n")
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.tumor)
        # read interactions.csv in the tumor directory
        interactions = pd.read_csv(os.path.join(tumor_dir, "interactions.csv"))
        # make dict with all genes in each module
        modules_t = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}
    else:
        print("===== Normal =====\n") 
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.normal)
        # read interactions.csv in the normal directory
        interactions = pd.read_csv(os.path.join(normal_dir, "interactions.csv"))
        # make dict with all genes in each module
        modules_n = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}

    genes = WGCNA.datExpr.var

    # Extract genes for complex A and B for each row
    complex_a_genes = interactions.apply(lambda row: list(set([row[f'interactor{i}'] for i in range(1, row['num_interactors_a'] + 1) if pd.notna(row[f'interactor{i}'])])), axis=1)
    complex_b_genes = interactions.apply(lambda row: list(set([row[f'interactor{i}'] for i in range(row['num_interactors_a'] + 1, row['num_interactors_a'] + row['num_interactors_b'] + 1) if pd.notna(row[f'interactor{i}'])])), axis=1)
    interactions_all_genes = complex_a_genes + complex_b_genes

    # make all_interacting_genes a list containing all unique genes in the interactions
    all_interacting_genes = list(set().union(*interactions_all_genes))
    print("Number of genes in CPDB: ", len(all_interacting_genes))
    
    # ----- ROC all gene pairs -----
    
    print("Generating ROC curves")
    
    # subset to all_interacting_genes
    tom = WGCNA.TOM
    #print("Subsetting to all interacting genes")
    #tom = WGCNA.TOM.loc[all_interacting_genes, all_interacting_genes]
    print("TOM has shape: ", tom.shape)
    print(tom.head())

    # print number of nan values
    print("Number of nans in TOM: ", tom.isna().sum().sum())

    # flatten TOM matrix, remove diagonal and duplicated values
    print("Flattening TOM matrix")
    all_pairs = pd.DataFrame(
        tom.where(
            np.tri(
                tom.shape[0],
                dtype=bool,
                k=-1
            ),
            np.nan
        ).stack(dropna=True), columns=["TOM"]
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
    # print nan rows
    print(all_pairs[all_pairs.isna().any(axis=1)])


    if condition == "tumor":
        name = WGCNA.name
    else:
        name = None

    # Diff complex interactions
    generate_roc_curve(
        data=all_pairs,
        target_col="interaction",
        feature_col="TOM",
        condition=condition,
        comparison_file=comparison_file,
        ax=roc_ax,
        name=name
    )
    print("Done: " + condition + " ROC curve")
    
    # ------ Test VS All -------
    print("Performing Mann-Whitney U test")
    
    # Perform test: do LR pairs have higher TOM than non-LR pairs?
    U_all, p_all = mannwhitneyu(
        all_pairs.loc[all_pairs['interaction'] == True, "TOM"],
        all_pairs.loc[all_pairs['interaction'] == False, "TOM"],
        alternative="greater"
    )
    print("Writing result to comparison.txt")
    
    # create it if it does not exist, overwrite it if it does
    with open(comparison_file, "a") as f:
        f.write(condition + "_mannwhitneyu_U_all: " + str(U_all) + "\n")
        f.write(condition + "_mannwhitneyu_p_all: " + str(p_all) + "\n")
        
    print("Done: " + condition + " test") 

# Save the roc figure
roc_fig.savefig(os.path.join(tumor_fig_dir, "roc.pdf"))

# --------- Jaccard Similarity---------
print("Calculating Jaccard similarity")

df_n = pd.read_csv(os.path.join(normal_dir, "interactions.csv"))
df_t = pd.read_csv(os.path.join(tumor_dir, "interactions.csv"))

interactions_n = df_n['interaction']
interactions_t = df_t['interaction']

same_module_n = df_n.loc[df_n["same_module"], "interaction"]
same_module_t = df_t.loc[df_t["same_module"], "interaction"]

def jaccard_similarity(set1, set2):
    intersection_size = len(set1.intersection(set2))
    union_size = len(set1.union(set2))
    similarity = intersection_size / union_size
    return similarity

# Calculate Jaccard similarity
js_all = jaccard_similarity(set(interactions_n), set(interactions_t))
js_same_module = jaccard_similarity(set(same_module_n), set(same_module_t))

# Print the result
print("Jaccard Similarity:")
print(f"All LR pairs: {js_all}")
print(f"Same module: {js_same_module}")

# Add the result to the comparison.txt file in the tumor directory 
with open(comparison_file, "a") as f:
    f.write("jaccard_all_interactions: " + str(js_all) + "\n")
    f.write("jaccard_same_module_interactions: " + str(js_same_module) + "\n")

# ----- Module Similarity -----
print("Calculating module similarity")

# For each pair of one Normal and one Tumor module, calculate the Jaccard similarity
# Save in a dataframe with the module names as index and columns
module_similarity = pd.DataFrame(index=modules_n.keys(), columns=modules_t.keys())

for n in modules_n.keys():
    for t in modules_t.keys():
        js = jaccard_similarity(set(modules_n[n]), set(modules_t[t]))
        module_similarity.loc[n, t] = js

# Prepend N to the row names 
module_similarity.index = "N" + module_similarity.index.astype(str)
# Prepend T to the column names
module_similarity.columns = "T" + module_similarity.columns.astype(str)

# Save the result in the tumor directory
print("Saving module similarity to: " + os.path.join(tumor_dir, "module_similarity.csv"))
module_similarity.to_csv(os.path.join(tumor_dir, "module_similarity.csv"))

print("Done: all")
