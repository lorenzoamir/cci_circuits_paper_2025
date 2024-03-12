import numpy as np
import sys
import scanpy as sc
import PyWGCNA
import pandas as pd
import os
import argparse
from pydeseq2.dds import DeseqDataSet

# read filename from command line
#filename = os.environ.get("WGCNA_FILE_PATH")

parser = argparse.ArgumentParser(description='Run WGCNA pipeline')

parser.add_argument('-i', '--input', type=str, help='path to input adata file')
parser.add_argument('-n', '--ncpus', type=int, help='number of cpus to use')
parser.add_argument('-r', '--lr_resource', type=str, help='path to ligand-receptor resource')

args = parser.parse_args()

filename = args.input
adata = sc.read_h5ad(filename)
n_cpus = args.ncpus

# category is the last part of the filename after the last /
category = filename.split("/")[-1].split(".")[0]
# output_path is the path to the category folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("category", category)
print("output_path", output_path)
print("filename", filename)

print("Loading LR resource: {}".format(args.lr_resource))
resource = pd.read_csv(args.lr_resource)

# --------- Check corresponding normal tissue -----------
print("Checking corresponding normal tissue")
normal_tissue_files = None

if "/tumor/" in filename:
    print("Checking corresponding normal tissue")
    # Get the corresponding normal tissue, strip everything after /tumor/
    normal_tissue = filename.split("/tumor/")[0] + "/normal/"
    # Search for .h5ad files in the normal tissue folder
    # and save their absolute path
    normal_tissue_files = [os.path.join(normal_tissue, f) for f in os.listdir(normal_tissue) if f.endswith(".h5ad")]
    # replace .h5ad with .p
    normal_tissue_files = [f.replace(".h5ad", ".p") for f in normal_tissue_files]
    # add "WGCNA_" prefix
    normal_tissue_files = [f.replace("/normal/", "/normal/WGCNA_") for f in normal_tissue_files]
    # Check if there is only one corresponding normal tissue
    if len(normal_tissue_files) == 1:
        print("Found corresponding normal tissue: {}".format(normal_tissue_files[0]))
        # Writing normal tissue location in same folder as tumor tissue
        with open(os.path.join(output_path, "corresponding_normal_file.txt"), "w") as f:
            f.write(normal_tissue_files[0])
    elif len(normal_tissue_files) == 0:
        print("No corresponding normal tissue found")
        normal_tissue_files = None
    else:
        print("Multiple corresponding normal tissues found: {}".format(normal_tissue_files))
        normal_tissue_files = None

print("Done: Checking corresponding normal tissue")

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
dds = DeseqDataSet(
    counts=pd.DataFrame(data=adata.X.astype(int), index=adata.obs_names, columns=adata.var_names),
    metadata=adata.obs,
    design_factors=design_factors,  # compare samples based on the "condition"
    refit_cooks=True,
    #n_cpus = n_cpus
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

print("Done: DeSeq2")

# --------- WGCNA ------------

print("Starting WGCNA")

# make expression dataframe 
geneExp = pd.DataFrame(data=adata.X,  index=adata.obs_names, columns=adata.var_names)

WGCNA = PyWGCNA.WGCNA(
    name="WGCNA_" + category, 
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

n_modules = WGCNA.datExpr.var["moduleLabels"].nunique()
module_sizes = WGCNA.datExpr.var["moduleLabels"].value_counts().to_dict()

print("Done: WGCNA")

# Get all module names


# ------ TOM & adjacency ------

print("Saving TOM and adjacency")

# save TOM matrix as csv
WGCNA.TOM.to_csv(os.path.join(output_path, "tom.csv.gz"))
# save adjacency matrix as csv
WGCNA.adjacency.to_csv(os.path.join(output_path, "adjacency.csv.gz"))

print("Done: TOM and adjacency")

# -------- Same Module ---------

print("Checking which gene pairs are in the same module")

same_module = WGCNA.datExpr.var['moduleLabels'].values[:, None] == WGCNA.datExpr.var['moduleLabels'].values
same_module = pd.DataFrame(data=same_module, index=WGCNA.datExpr.var.index, columns=WGCNA.datExpr.var.index)

# save same_module as csv
same_module.to_csv(os.path.join(output_path, "same_module.csv.gz"))

# Get idxs corresponding to upper triangle (diagonal excluded)
idxs_x, idxs_y = np.triu_indices(same_module.shape[0], 1)

same_module = same_module.values[idxs_x, idxs_y]

total_pairs = same_module.shape[0]
n_same_module = same_module.sum()

# --------- LR Pairs ---------

genes = WGCNA.datExpr.var

print("Annotating genes")

# Drop interactions in which the ligand or the receptor is missing from the WGCNA data
resource = resource[resource["ligand"].isin(genes.index) & resource["receptor"].isin(genes.index)]
print("Only keeping interactions with both ligand and receptor in WGCNA data: {} interactions".format(resource.shape[0]))

# For each ligand-receptor check whether the ligand and receptor are in the same WGCNA module
resource["same_module"] = False
resource["L_module"] = None
resource["R_module"] = None
for i, row in resource.iterrows():
    ligand_module = genes.loc[row["ligand"], "moduleLabels"]
    receptor_module = genes.loc[row["receptor"], "moduleLabels"]
    # Add L_module and R_module columns containing moduleLabels
    resource.loc[i, "L_module"] = ligand_module
    resource.loc[i, "R_module"] = receptor_module
    if ligand_module == receptor_module:
        resource.loc[i, "same_module"] = True

# Add TOM to LR interactions
print("Adding TOM to LR interactions")
tom = WGCNA.TOM
tom["gene"] = tom.index
merged_df = pd.merge(resource, tom, left_on='ligand', right_on='gene').drop('gene', axis=1)
resource['TOM'] = merged_df.apply(lambda row: tom.loc[row['receptor'], row['ligand']], axis=1)

# Add adjacency to LR interactions
print("Adding adjacency to LR interactions")
adj = WGCNA.adjacency
adj["gene"] = adj.index
merged_df = pd.merge(resource, adj, left_on='ligand', right_on='gene').drop('gene', axis=1)
resource['adjacency'] = merged_df.apply(lambda row: adj.loc[row['receptor'], row['ligand']], axis=1)

print("Showing the first 10 interactions")
print(resource.head(10))

# Save the interactions in the same path as the WGCNA object

print("Saving interactions to {}".format(os.path.join(output_path, "LR_interactions.csv")))
resource.to_csv(os.path.join(output_path, "LR_interactions.csv"), index=False)

print("Done: LR Pairs")

# Z-score number of LR pairs in the same module using a binomial
# with n = total_LR_pairs and p = n_same_module / total_pairs

print("Z-scoring number of LR pairs in the same module")

n = resource.shape[0]
p = n_same_module / total_pairs

n_lr_same_module = resource["same_module"].sum()

z_score = (n_lr_same_module - n * p) / np.sqrt(n * p * (1 - p))

# ----- General info ------

print("Saving general info")
with open(os.path.join(output_path, "general_info.txt"), "w") as f:
    f.write("name: {}\n".format("WGCNA_" + category))
    f.write("WGCNA_file: {}\n".format(os.path.join(output_path, "WGCNA_" + category + ".p")))
    f.write("LR_resource: {}\n".format(args.lr_resource))
    if normal_tissue_files:
        if len(normal_tissue_files) == 1:
            f.write("corresponding_normal_tissue: {}\n".format(normal_tissue_files[0]))

    f.write("n_modules: {}\n".format(n_modules))
    f.write("module_sizes: {}\n".format(module_sizes))

    f.write("n_genes: {}\n".format(genes.shape[0]))
    f.write("n_samples: {}\n".format(WGCNA.datExpr.shape[0]))

    f.write("total_gene_pairs: {}\n".format(total_pairs))
    f.write("gene_pairs_same_module: {}\n".format(n_same_module))

    f.write("total_LR_pairs: {}\n".format(resource.shape[0]))
    f.write("LR_pairs_same_module: {}\n".format(n_lr_same_module))
    f.write("LR_pairs_Z_score: {}\n".format(z_score))

print("Done: wgcna.py")

# ----- Done ------

print("Done: WGCNA pipeline")
