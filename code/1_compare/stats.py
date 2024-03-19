import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import PyWGCNA
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import ranksums
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

with open(info_file, "r") as f:
    for line in f:
        if "LR_resource" in line:
            resource_path = line.split(": ")[1].strip()
            break

print("Using resource file: " + resource_path)
resource = pd.read_csv(resource_path)
lr_genes = set(resource["ligand"]).union(resource["receptor"])

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
roc_fig_vs_all, roc_ax_vs_all = plt.subplots(1, 1, figsize=(10,10))
roc_fig_vs_lr, roc_ax_vs_lr = plt.subplots(1, 1, figsize=(10,10))

for condition in ["normal", "tumor"]:

    if condition == "tumor":
        print("===== Tumor =====\n")
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.tumor)
        # make dict with all genes in each module
        modules_t = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}
    else:
        print("===== Normal =====\n") 
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.normal)
        # make dict with all genes in each module
        modules_n = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}

    # ----------- Heatmap ------------
    print("Generating heatmap")

    ordered = WGCNA.datExpr.var.sort_values(by="moduleLabels").index
    tom_all = WGCNA.TOM.loc[ordered, ordered]
    
    lr_ordered = ordered[ordered.isin(lr_genes)]
    tom_lr = WGCNA.TOM.loc[lr_ordered, lr_ordered]
    
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    
    fig.suptitle(WGCNA.name, fontsize=fs)
    
    vmax = np.percentile(tom_lr, 99)
    vmin = 0
    cax = ax[0].imshow(tom_all, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    ax[0].set_title("All genes", fontsize=fs)
    cax = ax[1].imshow(tom_lr, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    ax[1].set_title("Ligands & receptors", fontsize=fs)
    
    plt.tight_layout()

    # save figure to the figure folder of the tumor or normal WGCNA

    if condition == "tumor":
        # get the directory of the tumor WGCNA
        output_path = os.path.join(tumor_fig_dir, "heatmap_all_genes_vs_lr.pdf")
    else:
        # get the directory of the normal WGCNA
        output_path = os.path.join(normal_fig_dir, "heatmap_all_genes_vs_lr.pdf")

    plt.savefig(output_path)

    print("Done: " + condition + " heatmap")
    
    # ----- ROC all gene pairs -----
    
    print("Generating ROC curve: all gene pairs vs LR interactions")

    # flatten TOM matrix, remove diagonal and duplicated values
    all_gene_pairs = pd.DataFrame(
        WGCNA.TOM.where(
            np.tri(
                WGCNA.TOM.shape[0],
                dtype=bool,
                k=-1
            ),
            np.nan
        ).stack(dropna=True), columns=["TOM"]
    )
    
    #all_gene_pairs = pd.DataFrame(
    #    WGCNA.TOM.where(
    #        np.tri(
    #            WGCNA.TOM.shape[0], dtype=bool, k=-1
    #        ),
    #        np.nan
    #    ).stack(dropna=True).reset_index(name="TOM")
    #)

    allgenes = set(all_gene_pairs.index.get_level_values(0))

    all_gene_pairs["LR Pair"] = False

    print("Total number of gene pairs: ", len(all_gene_pairs))
    
    # Identify rows where both ligand and receptor are present in allgenes
    valid_rows = resource[(resource['ligand'].isin(allgenes)) & (resource['receptor'].isin(allgenes))]

    # Set "LR Pair" to True for LR interaction pairs
    # There's no risk of duplicates because all_gene_pairs does not contain duplicates

    zip_lr = pd.Series(zip(valid_rows['ligand'], valid_rows['receptor']))
    zip_lr = zip_lr[zip_lr.isin(all_gene_pairs.index)]
    all_gene_pairs.loc[zip_lr, "LR Pair"] = True

    zip_rl = pd.Series(zip(valid_rows['receptor'], valid_rows['ligand']))
    zip_rl = zip_rl[zip_rl.isin(all_gene_pairs.index)]
    all_gene_pairs.loc[zip_rl, "LR Pair"] = True
    
    # Number of LR interactions:
    print("Number of LR interactions: ", all_gene_pairs["LR Pair"].sum())

    # X is the feature (TOM) and y is the target variable (LR Pair)

    # ROC curve using all gene pairs
    fpr_vs_all, tpr_vs_all, _ = roc_curve(all_gene_pairs['LR Pair'], all_gene_pairs['TOM'])
    roc_auc_vs_all = roc_auc_score(all_gene_pairs['LR Pair'], all_gene_pairs['TOM'])
    
    if condition == "tumor":
        with open(comparison_file, "a") as f:
            f.write("roc_auc_tumor_vs_all: " + str(roc_auc_vs_all) + "\n")
    elif condition == "normal":
        with open(comparison_file, "a") as f:
            f.write("roc_auc_normal_vs_all: " + str(roc_auc_vs_all) + "\n")

    # Plot ROC curve and include AUC in the legend
    if condition == "tumor":
        color = 'C1' 
        label = f'TCGA {WGCNA.name} (AUC = {roc_auc_vs_all:.2f})'

    elif condition == "normal":
        color = 'C0'
        label = f'GTEx {WGCNA.name} (AUC = {roc_auc_vs_all:.2f})'

    roc_ax_vs_all.plot([0, 1], [0, 1], color='k', linestyle='--', lw=2)
    roc_ax_vs_all.plot(fpr_vs_all, tpr_vs_all, color=color, lw=2, label=label)
    roc_ax_vs_all.set_xlabel('False Positive Rate')
    roc_ax_vs_all.set_ylabel('True Positive Rate')
    roc_ax_vs_all.set_title(WGCNA.name + " LR interactions vs all gene pairs", fontsize=fs)
    roc_ax_vs_all.legend()

    # save roc figure to the figure folder of the tumor WGCNA
    output_path = os.path.join(tumor_fig_dir, "roc_vs_all.pdf")
    roc_fig_vs_all.savefig(output_path)

    print("Done: " + condition + " ROC curve vs all gene pairs")

    # ----- ROC LR -----

    print("Generating ROC curve: LR combinations")
    
    all_lr_combinations = all_gene_pairs.loc[
        all_gene_pairs.index.get_level_values(0).isin(lr_genes) & \
        all_gene_pairs.index.get_level_values(1).isin(lr_genes) \
    ]

    fpr_vs_lr, tpr_vs_lr, _ = roc_curve(
        all_lr_combinations['LR Pair'],
        all_lr_combinations['TOM']
    )

    roc_auc_vs_lr = roc_auc_score(
        all_lr_combinations['LR Pair'],
        all_lr_combinations['TOM']
    )

    if condition == "tumor":
        with open(comparison_file, "a") as f:
            f.write("roc_auc_tumor_vs_lr: " + str(roc_auc_vs_lr) + "\n")
    elif condition == "normal":
        with open(comparison_file, "a") as f:
            f.write("roc_auc_normal_vs_lr: " + str(roc_auc_vs_lr) + "\n")

    # Plot ROC curve and include AUC in the legend
    if condition == "tumor":
        color = 'C1' 
        label = f'TCGA {WGCNA.name} (AUC = {roc_auc_vs_lr:.2f})'

    elif condition == "normal":
        color = 'C0'
        label = f'GTEx {WGCNA.name} (AUC = {roc_auc_vs_lr:.2f})'

    roc_ax_vs_lr.plot([0, 1], [0, 1], color='k', linestyle='--', lw=2)
    roc_ax_vs_lr.plot(fpr_vs_lr, tpr_vs_lr, color=color, lw=2, label=label)
    roc_ax_vs_lr.set_xlabel('False Positive Rate')
    roc_ax_vs_lr.set_ylabel('True Positive Rate')
    roc_ax_vs_lr.set_title(WGCNA.name + " LR combinations", fontsize=fs)
    roc_ax_vs_lr.legend()
    
    # save roc figure to the figure folder of the tumor WGCNA
    output_path = os.path.join(tumor_fig_dir, "roc_vs_lr.pdf")
    roc_fig_vs_lr.savefig(output_path)

    print("Done: " + condition + " ROC curve vs LR pairs")
    
    # ------ Rank Sum Test VS All -------
    print("Performing rank sum test vs all gene pairs")

    all_gene_pairs.loc[all_gene_pairs["LR Pair"] == False, "TOM"]
    lr_tom = all_gene_pairs.loc[all_gene_pairs["LR Pair"] == True, "TOM"]
    
    # Perform rank sum test: do LR pairs have higher TOM than non-LR pairs?
    U_vs_all, p_vs_all = ranksums(
        all_gene_pairs.loc[all_gene_pairs["LR Pair"] == True, "TOM"],
        all_gene_pairs.loc[all_gene_pairs["LR Pair"] == False, "TOM"],
        alternative="greater"
    )

    print("Writing result to comparison.txt")
    
    # create it if it does not exist, overwrite it if it does
    with open(comparison_file, "a") as f:
        f.write(condition + "_ranksums_U_vs_all: " + str(U_vs_all) + "\n")
        f.write(condition + "_ranksums_p_vs_all: " + str(p_vs_all) + "\n")
    print("Done: " + condition + " rank sum test vs all gene pairs") 

    # ------ Rank Sum Test VS LR -------

    print("Performing rank sum test vs LR pairs")
    # Rank sum ov LR interactions vs all ligand and receptors
    U_vs_lr, p_vs_lr = ranksums(
        all_lr_combinations.loc[all_lr_combinations["LR Pair"] == True, "TOM"],
        all_lr_combinations.loc[all_lr_combinations["LR Pair"] == False, "TOM"],
        alternative="greater"
    )
    
    print("Writing result to comparison.txt")
    
    with open(comparison_file, "a") as f:
        f.write(condition + "_ranksums_U_vs_lr: " + str(U_vs_lr) + "\n")
        f.write(condition + "_ranksums_p_vs_lr: " + str(p_vs_lr) + "\n")

    print("Done: " + condition + " rank sum test vs LR combinations")

    print()

# --------- Jaccard Similarity---------
print("Calculating Jaccard similarity")

LR_pairs_n = pd.read_csv(os.path.join(normal_dir, "LR_interactions.csv"))
LR_pairs_t = pd.read_csv(os.path.join(tumor_dir, "LR_interactions.csv"))

all_lr_pairs_n = set(LR_pairs_n["ligand"] + "-" + LR_pairs_n["receptor"])
all_lr_pairs_t = set(LR_pairs_t["ligand"] + "-" + LR_pairs_t["receptor"])

subset_n = LR_pairs_n.loc[LR_pairs_n["same_module"]]
subset_t = LR_pairs_t.loc[LR_pairs_t["same_module"]]

same_module_lr_pairs_n = set(subset_n["ligand"] + "-" + subset_n["receptor"])
same_module_lr_pairs_t = set(subset_t["ligand"] + "-" + subset_t["receptor"])

def jaccard_similarity(set1, set2):
    intersection_size = len(set1.intersection(set2))
    union_size = len(set1.union(set2))
    similarity = intersection_size / union_size
    return similarity

# Calculate Jaccard similarity
js_all = jaccard_similarity(all_lr_pairs_n, all_lr_pairs_t)
js_same_module = jaccard_similarity(same_module_lr_pairs_n, same_module_lr_pairs_t)

# Print the result
print("Jaccard Similarity:")
print(f"All LR pairs: {js_all}")
print(f"Same module: {js_same_module}")

# Add the result to the comparison.txt file in the tumor directory 
with open(comparison_file, "a") as f:
    f.write("jaccard_all_lr_pairs: " + str(js_all) + "\n")
    f.write("jaccard_same_module: " + str(js_same_module) + "\n")

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
