import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
#from combat import pycombat
from pycombat import Combat
import os
import re
import sys

seed = 42
BATCH_CORRECTION = True

# TIGER Data
data = pd.read_csv('/projects/bioinformatics/DB/tiger_immunotherapy/full_merged_dataset.csv', index_col=0)

# ERR2498029 is a clear outlier in the PCA plot
data = data.loc[data.index != 'ERR2498029']

# Convert to TPMs to match other data sources
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# http://arxiv.org/abs/1104.3889
# tpms = 1e6 * (data.iloc[:, 16:].T / data.iloc[:, 16:].sum(1)).T
# data.iloc[:, 16:] = tpms

# Match TCGA data log2(fpkm+0.001)
data.iloc[:, 16:] = np.log2(data.iloc[:, 16:] + 0.001)

# Rename to match gene names from ccc list
rename = rename = {
    # tiger  ->  ccc
    'ACPP': 'ACP3',
    'CECR1': 'ADA2',
    #'CMKLR2': 'GPR1', # We have both of them in the ccc list, so either name should be fine
    'CD8BP': 'CD8B2',
    'EPRS': 'EPRS1', # could also be QARS
    'HLA-DRB3': 'HLA-DPB1',
    'FAM213B': 'PRXL2B',
    'FAM19A1': 'TAFA1',
    'FAM19A4': 'TAFA4',
    'FAM19A5': 'TAFA5',
    'C10orf54': 'VSIR',
}
data = data.rename(columns=rename)

data = data.rename(columns={'Therapy': 'therapy_type'})

# Only keep RNA-Seq
data = data[data.seq_type.isin(['RNA-seq', 'RNA-Seq'])]
clinical_cols = list(data.columns[:16])

tumor_type_to_tissue = {
    'Melanoma' : 'Skin',
    'ccRCC' : 'Kidney',
    'Metastatic_gastric_cancer' : 'Stomach',
    'NSCLC' : 'Lung',
    'GBM' : 'Brain',
}

data['tissue'] = data['tumor_type'].map(tumor_type_to_tissue)
data['batch'] = 1
clinical_cols += ['tissue', 'batch']

mapping = {
    'PRE': 'PRE',
    'ON': 'ON',
    'POST': 'POST',
    'POST1': 'POST',
    'POST2': 'POST',
    'POST3': 'POST',
    'POST4': 'POST',
    'POST5': 'POST',
}

# Map values, the ones not in mapping will be NaN
data['treatment_when'] = data['Treatment'].map(mapping)
clinical_cols += ['treatment_when']

# TCGA Data

harmo = pd.read_csv('/home/lnemati/tcga_drug_harmonized.tsv', usecols=['measure_of_response', 'therapy_type', 'drug_name', 'tumor', 'patient_barcode'], sep='\t')
harmo = harmo.dropna(subset='measure_of_response')

harmo['response_NR'] = np.where(
    harmo.measure_of_response.isin(['complete response', 'partial response']),
    'R',
    'N'
)

immunotherapy_with_response = harmo.query('therapy_type == "immunotherapy"')['patient_barcode'].unique()

harmo['therapy_type'] = np.nan
harmo.loc[harmo['drug_name'] == "IPILIMUMAB", 'therapy_type'] = 'anti-CTLA-4'

ctla = harmo.dropna(subset='therapy_type').patient_barcode.unique()

# TPMs
#surv = pd.read_csv('/home/lnemati/pathway_crosstalk/data/survival_data/tissues/pan_cancer.csv', index_col=0)

# FPKMs
surv = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/tcga_log_fpkms.csv', index_col=0)

surv['patient_name'] = [x.rsplit('-', 1)[0] for x in surv.index]
surv['patient_name'] = surv['patient_name'].str.lower()

surv['dataset_id'] = [x.rsplit('-', 2)[0] for x in surv.index]
surv['dataset_id'] = surv['dataset_id'].str.lower()

surv['treatment_when'] = np.nan

patient_to_response = harmo.loc[harmo.patient_barcode.isin(immunotherapy_with_response), ['patient_barcode', 'response_NR']]
patient_to_response = patient_to_response.drop_duplicates()
patient_to_response = patient_to_response.set_index('patient_barcode')['response_NR']

surv['response_NR'] = surv['patient_name'].replace(patient_to_response)
surv.loc[~surv['response_NR'].isin(['R', 'N']), 'response_NR'] = np.nan

surv['therapy_type'] = np.nan
surv.loc[surv['patient_name'].isin(ctla), 'therapy_type'] = 'anti-CTLA-4'
surv['batch'] = 2

common_cols = surv.columns.intersection(data.columns)
data = data[common_cols]
surv = surv[common_cols]

print('Merging TIGER and TCGA data')
data = pd.concat([data, surv])

# IMVigor210 Data
im = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/imvigor210/merged_log_fpkm.csv', index_col=0)
im['batch'] = 3

common_cols = im.columns.intersection(data.columns)
data = data[common_cols]
im = im[common_cols]

print('Merging TIGER, TCGA and IMVigor210 data')
data = pd.concat([data, im])
print(data.shape)

data.loc[~data['response_NR'].isin(['R', 'N']), 'response_NR'] = np.nan

print(data.response_NR.value_counts(dropna=False))

# Discard patients with many missing values
#data = data[data.isna().sum(1) < 3000]

# Make sure patient names are unique
data.patient_name = data['dataset_id'].astype(str) + '-' + data.patient_name.astype(str)

clinical_cols = list(set(clinical_cols).intersection(common_cols)) 
genes = [col for col in data.columns if col not in clinical_cols]

## BATCH CORRECTION (OLD)
#
## Remove batches with too few patients
#small_batch = data.dataset_id.value_counts()[data.dataset_id.value_counts() <= 5].index
#data = data[~data.dataset_id.isin(small_batch)]
#batch = list(data.dataset_id)
#
#expr = data[genes].T.fillna(0)
#corrected = pycombat.pycombat(data=expr, batch=batch)
#data[genes] = corrected.T

# Remove patient missing response data
print(data.response_NR.value_counts(dropna=False))
data = data.dropna(subset=['response_NR'])
print(data.response_NR.value_counts(dropna=False))

# Remove post-treatment samples
data = data[data['treatment_when'] != 'POST']

# Extract unique patients with known response
patients = data.loc[
    :,
    ['patient_name', 'tissue', 'response_NR']
].drop_duplicates()

# TRAIN TEST SPLIT
print('Splitting patients')
train_patients, test_patients = train_test_split(
    patients,
    test_size=0.2,
    random_state=seed
)

# Filter the original data based on the selected patients
test = data[data['patient_name'].isin(test_patients['patient_name'])]
train = data[~data['patient_name'].isin(test_patients['patient_name'])]

assert set(train['patient_name']).isdisjoint(set(test['patient_name'])), "Error: Overlapping patients in train and test!"

# BATCH CORRECTION (NEW)
if BATCH_CORRECTION:
    print('Batches:')
    print(train.batch.value_counts())
    print(test.batch.value_counts())

    def preprocess_data(data, genes):
        """Extract gene expression matrix and phenotype data."""
        X = data[genes].fillna(np.log2(0.001))
        obs = data.drop(columns=genes)
        return X.values, obs

    def encode_response(pheno):
        """Encode response variable as a binary matrix."""
        return pd.get_dummies(pheno.response_NR, drop_first=True).values

    def apply_combat(train_Y, train_b, train_X, test_Y, test_b, test_X, train_C=None, test_C=None):
        """Apply Combat batch effect correction."""
        combat = Combat()
        combat.fit(Y=train_Y, b=train_b, X=train_X, C=train_C)
        return combat.transform(Y=train_Y, b=train_b, X=train_X, C=train_C), combat.transform(Y=test_Y, b=test_b, X=test_X, C=test_C)

    # Preprocess train and test data
    train_Y, train_pheno = preprocess_data(train, genes)
    test_Y, test_pheno = preprocess_data(test, genes)

    # Encode response
    train_X = encode_response(train_pheno)
    test_X = encode_response(test_pheno)

    # Extract batch values
    train_b = train_pheno['batch'].values
    test_b = test_pheno['batch'].values

    # Apply Combat correction
    train_combat, test_combat = apply_combat(train_Y, train_b, train_X, test_Y, test_b, test_X)

    # Update dataframes
    train[genes] = train_combat
    test[genes] = test_combat

train.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/train.csv')
test.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/test.csv')

print('Train:', train.shape)
print('Test:', test.shape)

print('Done: split.py')
