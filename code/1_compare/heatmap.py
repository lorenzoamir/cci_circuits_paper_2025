import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import PyWGCNA
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Compare WGCNA results')

parser.add_argument('--tumor', type=str, help='Input tumor WGCNA')
parser.add_argument('--normal', type=str, help='Input normal WGCNA')

args = parser.parse_args()

print("normal file:", args.normal)
print("tumor file:", args.tumor)

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

# Function that generates and saves the ROC curve
def generate_roc_curve(

for condition in ["normal", "tumor"]:

    if condition == "tumor":
        print("===== Tumor =====\n")
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.tumor)
        # make dict with all genes in each module
        modules_t = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}
        # read interactions.csv in the tumor directory
        cpdb = pd.read_csv(os.path.join(tumor_dir, "interactions.csv"))

    else:
        print("===== Normal =====\n") 
        print("Reading WGCNA ojbects")
        WGCNA = PyWGCNA.readWGCNA(args.normal)
        # make dict with all genes in each module
        modules_n = {i: WGCNA.datExpr.var.loc[WGCNA.datExpr.var["moduleLabels"] == i].index for i in WGCNA.datExpr.var["moduleLabels"].unique()}
        # read interactions.csv in the normal directory
        cpdb = pd.read_csv(os.path.join(normal_dir, "interactions.csv"))

    # Extract genes for complex A and B for each row
    complex_a_genes = cpdb.apply(lambda row: [row[f'interactor{i}'] for i in range(1, row['num_interactors_a'] + 1) if pd.notna(row[f'interactor{i}'])], axis=1)
    complex_b_genes = cpdb.apply(lambda row: [row[f'interactor{i}'] for i in range(row['num_interactors_a'] + 1, row['num_interactors_a'] + row['num_interactors_b'] + 1) if pd.notna(row[f'interactor{i}'])], axis=1)
    interactions_all_genes = complex_a_genes + complex_b_genes

    all_cpdb_genes = set([gene for genes in interactions_all_genes for gene in genes])
    print("Number of genes in CPDB: ", len(all_cpdb_genes))

    genes = WGCNA.datExpr.var

    # Only keep interactions where all genes are in the WGCNA data
    print(f'Number of interactions before filtering: {len(cpdb)}')
    cpdb = cpdb[interactions_all_genes.apply(lambda x: all([gene in genes.index for gene in x]))] 
        print("Keeping interactions with all genes present: ", len(cpdb))

    # ----------- Heatmap ------------
    print("Generating heatmap")

    ordered = WGCNA.datExpr.var.sort_values(by="moduleLabels").index
    tom_all = WGCNA.TOM.loc[ordered, ordered]
    
    #lr_ordered = ordered[ordered.isin(lr_genes)]
    cpdb_ordered = oredered[ordered.isin(all_cpdb_genes)]
    tom_cpdb = WGCNA.TOM.loc[cpdb_ordered, cpdb_ordered]
    
    fig, ax = plt.subplots(1,2, figsize=(20,10))
    
    fig.suptitle(WGCNA.name, fontsize=fs)
    
    vmax = np.percentile(tom_cpdb, 99)
    vmin = 0
    cax = ax[0].imshow(tom_all, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    ax[0].set_title("All genes", fontsize=fs)
    cax = ax[1].imshow(tom_cpdb, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax)
    ax[1].set_title("CellPhoneDB genes", fontsize=fs)
    
    plt.tight_layout()

    # save figure to the figure folder of the tumor or normal WGCNA
    if condition == "tumor":
        # get the directory of the tumor WGCNA
        output_path = os.path.join(tumor_fig_dir, "heatmap_all_genes_vs_cpdb.pdf")
    else:
        # get the directory of the normal WGCNA
        output_path = os.path.join(normal_fig_dir, "heatmap_all_genes_vs_cpdb.pdf")

    plt.savefig(output_path)

    print("Done: " + condition + " heatmap")

print("Done: heatmap")
