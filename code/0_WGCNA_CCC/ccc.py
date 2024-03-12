import numpy as np
import cell2cell as c2c
import scanpy as sc
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Run CCC pipeline')

parser.add_argument('-i', '--input', type=str, help='path to input adata file')

args = parser.parse_args()

print("Reading adata object:", args.input)
filename = args.input
adata = sc.read_h5ad(filename)

# category is the last part of the filename after the last /
category = filename.split("/")[-1].split(".")[0]
# output_path is the path to the category folder
output_path = "/".join(filename.split("/")[:-1]) + "/"

print("category", category)
print("output_path", output_path)
print("filename", filename)

# --- Subsetting ---
print(f"adata object has {adata.n_obs} samples")

max_n_samples = 100

# If there are more than max_n_samples samples, subset to 100 random samples
if adata.n_obs > max_n_samples: 
    print("Subsetting to max_n_samples random samples")
    adata = adata[np.random.choice(adata.obs_names, max_n_samples, replace=False), :]

# --- Cell-Cell Communication (cell2cell) ---

print("Getting TPMs from adata object")
rna_tpm = pd.DataFrame(adata.layers["TPM"].T, index=adata.var_names, columns=adata.obs_names)
metadata = adata.obs
metadata["sample_id"] = metadata.index

# read generic_info.txt file to get the path of LR_resource file
info_dir = filename.split("/")[:-1]
info_file = os.path.join(info_dir, "generic_info.txt")
print("Trying to read", info_file)
with open(info_file, "r") as f:
    lines = f.readlines()
    for line in lines:
        if "LR_resource" in line:
            # path is after "LR_resource: "
            resource = line.split(": ")[1].strip()
            print("Found LR_resource:", resource)

print("Loading LR resource")
resource = pd.read_csv(resource)

interactions = c2c.analysis.BulkInteractions(
    rnaseq_data=rna_tpm,
    ppi_data=resource,
    metadata=metadata,
    interaction_columns=("ligand", "receptor"),
    complex_sep=('_'),
    communication_score="expression_thresholding",
    expression_threshold=10, # TPMs
    cci_score="bray_curtis",
    cci_type="undirected",
    sample_col="sample_id",
    group_col="",
    verbose=False,
)

print("Computing CCC scores")
interactions.compute_pairwise_communication_scores(verbose=False)
#interactions.compute_pairwise_cci_scores(verbose=False)

print("Computing number and fraction of interactions")
# LR pairs x sample pairs
ccc = interactions.interaction_space.interaction_elements["communication_matrix"]
# Number of sample pairs in which the LR pair is expressed for each LR pair
n_interactions = ccc.sum(axis=1)
# Fraction of sample pairs in which the LR pair is expressed for each LR pair
frac_interactions = ccc.mean(axis=1)

# Save the LR pairs with their communication scores
# index is a tuple of ligand and receptor, create two columns for them
print("Creating dataframe with CCC scores")
ccc_scores = pd.DataFrame(n_interactions, columns=["n_interactions"])
ccc_scores["ligand"] = [i[0] for i in ccc_scores.index]
ccc_scores["receptor"] = [i[1] for i in ccc_scores.index]
ccc_scores["frac_interactions"] = frac_interactions

print(f"Merging scores to {os.path.join(output_path, 'LR_interactions.csv')}")
# read lr_pairs
lr_df = pd.read_csv(os.path.join(output_path, "LR_interactions.csv"))
# if lr_df already has the columns, remove them
lr_df = lr_df.drop(columns=["n_interactions", "frac_interactions"], errors="ignore")
# merge with ccc_scores
lr_df = lr_df.merge(ccc_scores, on=["ligand", "receptor"], how="left")

# !!! Warning: the merge procedure can introduce nans as some
# genes are both ligands and receptors!

print("Warning: the merge procedure can introduce nans!")

# save to file
lr_df.to_csv(os.path.join(output_path, "LR_interactions.csv"), index=False)

# ----- Done ------

print("Done: CCC")
print()
