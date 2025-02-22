import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from combat import pycombat
import os
import re
import sys

seed = 42

# TIGER Data
data = pd.read_csv('/projects/bioinformatics/DB/tiger_immunotherapy/full_merged_dataset.csv', index_col=0)

# ERR2498029 is a clear outlier in the PCA plot
data = data.loc[data.index != 'ERR2498029']

# Convert to TPMs to match other data sources
# https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
# http://arxiv.org/abs/1104.3889

tpms = 1e6 * (data.iloc[:, 16:].T / data.iloc[:, 16:].sum(1)).T
data.iloc[:, 16:] = tpms

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
data['response_known'] = np.where(data.response_NR == 'UNK', 'unknown', 'known')

# Divide treatment time into pre and all other options
treatment_when_nans = data.Treatment.isna()
data['treatment_when'] = np.where(data.Treatment == 'PRE', 'PRE', 'NOT_PRE')
data['treatment_when'] = np.where(treatment_when_nans, np.nan, data['treatment_when'])

tumor_type_to_tissue = {
    'Melanoma' : 'Skin',
    'ccRCC' : 'Kidney',
    'Metastatic_gastric_cancer' : 'Stomach',
    'NSCLC' : 'Lung',
    'GBM' : 'Brain',
}

data['tissue'] = data['tumor_type'].map(tumor_type_to_tissue)

clinical_cols = list(data.columns[:16])
clinical_cols += ['response_known', 'treatment_when', 'tissue']

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

surv = pd.read_csv('/home/lnemati/pathway_crosstalk/data/survival_data/tissues/pan_cancer.csv', index_col=0)

surv['patient_name'] = [x.rsplit('-', 1)[0] for x in surv.index]
surv['patient_name'] = surv['patient_name'].str.lower()

surv['dataset_id'] = [x.rsplit('-', 2)[0] for x in surv.index]
surv['dataset_id'] = surv['dataset_id'].str.lower()

surv['treatment_when'] = np.nan

patient_to_response = harmo.loc[harmo.patient_barcode.isin(immunotherapy_with_response), ['patient_barcode', 'response_NR']]
patient_to_response = patient_to_response.drop_duplicates()
patient_to_response = patient_to_response.set_index('patient_barcode')['response_NR']

surv['response_NR'] = surv['patient_name'].replace(patient_to_response)
surv['response_NR'] = np.where(~surv['response_NR'].isin(['N', 'R']), 'UNK', surv['response_NR'])

surv['response_known'] = np.where(surv['response_NR'] == 'UNK', 'unknown', 'known')

surv['therapy_type'] = np.nan
surv.loc[surv['patient_name'].isin(ctla), 'therapy_type'] = 'anti-CTLA-4'

common_cols = surv.columns.intersection(data.columns)
data = data[common_cols]
surv = surv[common_cols]

data = pd.concat([data, surv])

# Discard patients with many missing values
data = data[data.isna().sum(1) < 5000]

# Make sure patient names are unique
data.patient_name = data['dataset_id'].astype(str) + '-' + data.patient_name.astype(str)

clinical_cols = list(set(clinical_cols).intersection(common_cols)) 

# BATCH CORRECTION

# Remove batches with too few patients
small_batch = data.dataset_id.value_counts()[data.dataset_id.value_counts() <= 5].index
data = data[~data.dataset_id.isin(small_batch)]

genes = [col for col in data.columns if col not in clinical_cols]
batch = list(data.dataset_id)

expr = data[genes].T.fillna(0)
corrected = pycombat.pycombat(expr, batch)
data[genes] = corrected.T

# TRAIN TEST SPLIT

for col in clinical_cols:
    data[col] = data[col].astype("category").cat.remove_unused_categories()

# Extract unique patients with known response
strat_col = 'response_NR'
patients = data[data['response_known'] == "known"][['patient_name', 'tissue', 'response_NR']].drop_duplicates()

# Split unique patients into train and test sets while stratifying by tissue
print('Splitting patients')
train_patients, test_patients = train_test_split(
    patients,
    test_size=0.1,
    #stratify=patients[strat_col],
    random_state=seed
)

# Filter the original data based on the selected patients
test = data[data['patient_name'].isin(test_patients['patient_name'])]
test = test[test['response_known'] == 'known']
train = data[~data['patient_name'].isin(test_patients['patient_name'])]
train = train[train['response_known'] == 'known']

assert set(train['patient_name']).isdisjoint(set(test['patient_name'])), "Error: Overlapping patients in train and test!"

train.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/train.csv')
test.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/test.csv')

print('Train:', train.shape)
print(train.groupby('tissue')['response_NR'].value_counts())
print('Test:', test.shape)
print(test.groupby('tissue')['response_NR'].value_counts())

print('Done: split.py')
