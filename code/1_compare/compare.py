import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import PyWGCNA
from sklearn.metrics import roc_curve, roc_auc_score
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumor', type=str, help='Input tumor WGCNA')
parser.add_argument('--normal', type=str, help='Input normal WGCNA')

args = parser.parse_args()

print("normal file:", args.normal)
print("tumor file:", args.tumor)

resource_path = "/projects/bioinformatics/DB/CellCellCommunication/LianaResources/consensus.csv"
print("Using resource file: " + resource_path)
resource = pd.read_csv("/projects/bioinformatics/DB/CellCellCommunication/LianaResources/consensus.csv")
lr_genes = set(resource["ligand"]).union(resource["receptor"])

print("Reading WGCNA ojbects")
WGCNA_n = PyWGCNA.readWGCNA(args.normal)
WGCNA_t = PyWGCNA.readWGCNA(args.tumor)

# Plot heatmap of TOM for tumor and normal comparing all genes with ligand and receptors only
fs = 18
cmap='viridis'

# Create figure for the ROC curve
roc_fig, roc_ax = plt.subplots(1, 1, figsize=(10,10))

for condition in ["normal", "tumor"]:
    if condition == "tumor":
        WGCNA = WGCNA_t
    else:
        WGCNA = WGCNA_n
    
    ordered = WGCNA.datExpr.var.sort_values(by="moduleLabels").index
    tom_all = WGCNA.TOM.loc[ordered, ordered]
    
    lr_ordered = ordered[ordered.isin(lr_genes)]
    tom_lr = WGCNA.TOM.loc[lr_ordered, lr_ordered]
    
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    
    fig.suptitle(WGCNA.name, fontsize=fs)
    
    vmax = np.percentile(tom_lr, 99)
    vmin = 0
    cax = ax[0].imshow(tom_all, cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    ax[0].set_title("All genes", fontsize=fs)
    cax = ax[1].imshow(tom_lr, cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
    ax[1].set_title("Ligands & receptors", fontsize=fs)
    
    plt.tight_layout()

    # save figure to the figure folder of the tumor or normal WGCNA
    if condition == "tumor":
        # get the directory of the tumor WGCNA
        tumor_dir = os.path.join(os.path.dirname(args.tumor), "figures")
        output_path = os.path.join(tumor_dir, "heatmap_all_genes_vs_lr.pdf")
    else:
        # get the directory of the normal WGCNA
        normal_dir = os.path.join(os.path.dirname(args.normal), "figures")
        output_path = os.path.join(normal_dir, "heatmap_all_genes_vs_lr.pdf")
    plt.savefig(output_path)
    
    # ROC curve
    
    # flatten TOM matrix
    all_pairs = pd.DataFrame(WGCNA.TOM.stack(dropna=True, sort=True), columns=["TOM"])
    allgenes = set(all_pairs.index.get_level_values(0))
    all_pairs["LR Pair"] = False

    # Identify rows where both ligand and receptor are present in allgenes
    valid_rows = resource[(resource['ligand'].isin(allgenes)) & (resource['receptor'].isin(allgenes))]

    # Set "LR Pair" to True for the valid rows using loc and the index
    all_pairs.loc[list(zip(valid_rows['ligand'], valid_rows['receptor'])), "LR Pair"] = True

    # X is the feature (TOM) and y is the target variable (LR Pair)
    X = all_pairs['TOM']
    y = all_pairs['LR Pair']

    # Calculate ROC curve for different thresholds
    fpr, tpr, thresholds = roc_curve(y, X)

    # Calculate AUC (Area Under the Curve)
    roc_auc = roc_auc_score(y, X)

    # Plot ROC curve and include AUC in the legend
    if condition == "tumor":
        color = 'C1' 
        label = f'TCGA {WGCNA.name} (AUC = {roc_auc:.2f})'

    elif condition == "normal":
        color = 'C0'
        label = f'GTEx {WGCNA.name} (AUC = {roc_auc:.2f})'

    roc_ax.plot([0, 1], [0, 1], color='k', linestyle='--', lw=2)
    roc_ax.plot(fpr, tpr, color=color, lw=2, label=label)
    roc_ax.set_xlabel('False Positive Rate')
    roc_ax.set_ylabel('True Positive Rate')
    roc_ax.set_title(WGCNA.name, fontsize=fs)
    roc_ax.legend()

# save figure to the figure folder of the tumor WGCNA
output_path = os.path.join(tumor_dir, "roc_curve.pdf")
roc_fig.savefig(output_path)

print("Done")




