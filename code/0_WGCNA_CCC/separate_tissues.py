import pandas as pd
import os
import numpy as np
import scanpy as sc
import argparse

parser = argparse.ArgumentParser(description='Read data and separate it according to tissue and condition')

parser.add_argument('--inputfile', type=str, help='path to input adata file')
parser.add_argument('--outputdir', type=str, help='path to output adata file')

args = parser.parse_args()

print("Reading adata from {}".format(args.inputfile))
adata = sc.read_h5ad(args.inputfile)
print("adata has shape {}".format(adata.shape))

print("Subsetting samples")

rename = {
    "X_study": "study",
    "X_gender": "gender",
    "X_sample_type": "type",
    "X_primary_site": "tissue",
    "primary_disease_or_tissue": "condition"
}

adata.obs.rename(columns=rename, inplace=True)
adata.obs.drop(columns="detailed_category", inplace=True)

print("Column names after renaming: {}".format(adata.obs.columns))

# subset to primary tumor and healty patients
#adata = adata[~adata.obs.tissue.isna()]
adata = adata[adata.obs.type.isin(["Primary Tumor", "Normal Tissue"])]

# some Testicular Germ Cell Tumor are missing the gender
adata.obs.loc[adata.obs["condition"] == "Testicular Germ Cell Tumor", "gender"] = "Male"

#drop samples with missing metadata
#adata = adata[~(adata.obs.isna().sum(axis=1) > 0), :]

# drop samples with missing tissue or condition
adata = adata[~(adata.obs.tissue.isna() | adata.obs.condition.isna())]

print("Setting adata.X to raw counts")
# Set adata.X to raw counts
adata.X = np.matrix(adata.layers["raw_counts"].todense(), dtype=int) 

print("Calculating QC metrics")
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Filter samples based on QC metrics
print("Filtering samples based on QC metrics")
adata.obs["th_outlier"] = False

# n_genes_by_counts:
#
# min 20000
# max 40000

adata.obs["th_outlier"] = np.where(
    (adata.obs["n_genes_by_counts"] < 20000) | (adata.obs["n_genes_by_counts"] > 40000),
    True,
    adata.obs["th_outlier"]
)

print("n_genes_by_counts: {}".format(adata.obs.th_outlier.value_counts()))

# total_counts
#
# min 10000000 (1e7)
# max 120000000 (1.2e8)

adata.obs["th_outlier"] = np.where(
    (adata.obs["total_counts"] < 1e7) | (adata.obs["total_counts"] > 1.2e8),
    True,
    adata.obs["th_outlier"]
)

print("total_counts: {}".format(adata.obs.th_outlier.value_counts()))

# pct_counts_in_top_100_genes
# max 85

adata.obs["th_outlier"] = np.where(
    (adata.obs["pct_counts_in_top_100_genes"] > 90),
    True,
    adata.obs["th_outlier"]
)

print("pct_counts_in_top_100_genes: {}".format(adata.obs.th_outlier.value_counts()))

def is_outlier(adata, metric: str, nmads: int):
    """
    identify outliers based on a MAD threshold

    params:
        adata: AnnData object
        metric: adata.obs column to use as metric
        nmads: number of MADs to use as threshold
    """
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * (M - M.mean()).abs().mean()) | \
              (M > np.median(M) + nmads * (M - M.mean()).abs().mean())
    return outlier

def QC_filter(adata, nmad=3, batch_key="sample", subset=False, verbose=False):
    '''
    Filters low-quality samples in single-cell RNA-seq data based on QC metrics.

    :param adata: Anndata object representing the unfiltered data with QC metrics.
    :param nmad: Number of median absolute deviations (MADs) used for filtering "log1p_total_counts", "log1p_n_genes_by_counts", and "pct_counts_in_top_20_genes".
    :param batch_key: Key in the observation metadata specifying the batch/sample information.
    :param subset: If True, subsets the Anndata object to exclude the identified outliers.
    :param verbose: If True, prints additional information during the filtering process.
    '''

    adata.obs["qc_outlier"] = False

    for batch_id in adata.obs[batch_key].unique():
        # Subset to batch
        batch_indices = adata.obs[batch_key] == batch_id
        batch = adata[batch_indices].copy()

        if verbose:
            print(f"Batch: {batch_id}")

        if nmad > 0:
            if verbose:
                print(f"Thresholding at {nmad} MADs: log1p_total_counts, log1p_n_genes_by_counts, pct_counts_in_top_20_genes")
            qc_outliers = (
                is_outlier(batch, "log1p_total_counts", nmad)
                | is_outlier(batch, "log1p_n_genes_by_counts", nmad)
                | is_outlier(batch, "pct_counts_in_top_100_genes", nmad)
            )
            adata.obs.loc[batch_indices, "qc_outlier"] = qc_outliers

        if verbose:
            print(adata.obs.qc_outlier.value_counts())
            print()

    print(f"Total number of samples: {adata.n_obs}")
    print(f"Outliers: {sum(list(adata.obs['qc_outlier']))}")

    if subset:
        adata = adata[~adata.obs.outlier].copy()
        print(f"Number of samples after filtering: {adata.n_obs}")

    return adata
      
adata = QC_filter(adata, nmad=5, batch_key="condition", subset=False, verbose=False)

print("th_outlier: {}".format(adata.obs.th_outlier.value_counts()))
print("qc_outlier: {}".format(adata.obs.qc_outlier.value_counts()))

adata = adata[~(adata.obs["th_outlier"] | adata.obs["qc_outlier"])]

print("Filtered adata object has {} samples and {} genes".format(adata.n_obs, adata.n_vars))

print("Filtering and annotating genes")

print("Annotating protein coding genes")
# retrieve gene symbol of protein coding genes
annot = pd.read_csv("~/resources/biomart/ensembl_to_symbol_filtered.csv.gz", index_col="Gene stable ID")

# remove ensembl version
adata.var["ensembl"] = adata.var_names.str.replace(r"\..*","", regex=True)
adata.var_names = adata.var["ensembl"].values

# subset to genes in biomart
adata = adata[:, adata.var_names.isin(annot.index)]

# add gene symbols
adata.var["symbol"] = annot["Gene name"]

# remove genes with missing symbol
adata = adata[:, ~adata.var.symbol.isna()]

# set index to gene symbol
adata.var_names = adata.var.symbol.values

# only keep genes with counts > 15 in al least 75% of samples (suggested by WGCNA on RNAseq FAQ)
print("Filtering genes with counts < 15 in at least 75% of samples")
def filter_genes(adata, min_counts=15, fraction=0.75):
    """
    Filter genes based on count data in adata.X matrix.

    Parameters:
    - adata: AnnData object
    - min_counts: Minimum count threshold
    - max_fraction: Minimum fraction of samples where a gene should have counts above min_counts

    Returns:
    - Filtered AnnData object
    """

    # Check which combinations gene and sample are above min_counts
    is_above_thr = adata.X >= min_counts

    # For each gene calculate fraction of samples where gene has counts above min_counts
    fraction_per_gene = np.sum(is_above_thr, axis=0) / adata.shape[0]

    # Subset to only genes that have enough counts in a big enough fraction of samples
    return adata[:, fraction_per_gene > fraction]


# Divide according to tissue and condition

print("New data will be saved in {}".format(args.outputdir))
parent_dir = args.outputdir 

if not os.path.exists(parent_dir):
        os.makedirs(parent_dir)

## save a copy of the data without subsetting
#bdata = adata.copy()
#
## filtering genes
#bdata = filter_genes(bdata, min_counts=15, fraction=0.75)
#
## make All_Tissues directory
#tissue_path = os.path.join(parent_dir, "all_tissues")
#if not os.path.isdir(tissue_path):
#    os.mkdir(tissue_path)
#    print(f"Creating directory: {tissue_path}")
#
#filename = os.path.join(tissue_path, "all_tissues.h5ad")
#print(f"Writing file: {filename}")
#bdata.write(filename, compression="gzip")

print("columns are {}".format(adata.obs.columns))

print("Separating data according to tissue and condition")
print()

print("Tissues:")
print(adata.obs.tissue.value_counts())

print("Conditions:")
print(adata.obs.condition.value_counts())

def check_filter_write(adata, filename):
    """
    Filter genes, check if there are enough samples and genes left and write file 
    """
    # filtering genes
    if adata.shape[0] < 15:
        print(f"Less than 15 samples left after filtering. Skipping file: {filename}")
        return

    adata = filter_genes(adata, min_counts=15, fraction=0.75)
    if adata.shape[1] < 3000:
        print(f"Less than 3000 genes left after filtering. Skipping file: {filename}") 
        return

    print(f"Writing file: {filename}")
    adata.write(filename, compression="gzip")
    return

def clean_string(s):
    # substitute spaces and commas with underscores
    s = s.replace(" ", "_")
    s = s.replace(",", "_")
    # make lowercase
    s = s.lower()
    # remove non-alphanumeric characters except _ and &
    s = "".join([c for c in s if c.isalnum() or c in ["_", "&"]])
    # substitute __ with single _
    s = s.replace("__", "_")
    return s

for tissue in adata.obs.tissue.unique():

    print(tissue)

    # Path
    tissue_str = clean_string(tissue) 
    tissue_path = os.path.join(parent_dir, tissue_str)

    if not os.path.isdir(tissue_path):
        os.mkdir(tissue_path)
        print(f"Creating directory: {tissue_path}")

    bdata = adata[adata.obs.tissue== tissue]

    for condition in ["normal", "tumor"]:
        print(4*" " + condition)
        condition_path = os.path.join(tissue_path, condition)
        if not os.path.isdir(condition_path):
            os.mkdir(condition_path)
            print(f"Creating directory: {condition_path}")

        if condition == "normal":
            cdata = bdata[bdata.obs.type == "Normal Tissue"]
            filename = os.path.join(condition_path, tissue_str)
            filename += ".h5ad"
            print(f"Writing file: {filename}")
            check_filter_write(cdata, filename)

        elif condition == "tumor":
            cdata = bdata[bdata.obs.type == "Primary Tumor"]

            for tumor_type in cdata.obs.condition.unique():
                print(8*" " + tumor_type)
                tumor_type_str = clean_string(tumor_type)
                tumor_type_path = os.path.join(condition_path, tumor_type_str)
                if not os.path.isdir(tumor_type_path):
                    os.mkdir(tumor_type_path)
                    print(8*" "+f"Creating directory: {tumor_type_path}")
                ddata = cdata[cdata.obs.condition == tumor_type]
                filename = os.path.join(tumor_type_path, tumor_type_str)
                filename += ".h5ad"
                print(8*" "+f"Writing file: {filename}")
                check_filter_write(ddata, filename)
    print()

print("SeparateTissuesConditions.py finished")
