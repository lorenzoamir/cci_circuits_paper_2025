import pandas as pd
import numpy as np
from anndata import AnnData
from scanpy.pp import combat
import os
import re
import sys

seed = 42
BATCH_CORRECTION = True
np.random.seed(seed)

# TIGER Data
data = pd.read_csv('/projects/bioinformatics/DB/tiger_immunotherapy/full_merged_dataset.csv', index_col=0)

# ERR2498029 is a clear outlier in the PCA plot
#data = data.loc[data.index != 'ERR2498029']

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
clinical_cols += ['tissue']

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
surv = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/tcga/tcga_log_fpkms.csv', index_col=0)

surv['patient_name'] = [x.rsplit('-', 1)[0] for x in surv.index]
surv['patient_name'] = surv['patient_name'].str.lower()

#surv['dataset_id'] = [x.rsplit('-', 2)[0] for x in surv.index]
#surv['dataset_id'] = surv['dataset_id'].str.lower()
surv['dataset_id'] = 'tcga'

surv['treatment_when'] = np.nan

patient_to_response = harmo.loc[harmo.patient_barcode.isin(immunotherapy_with_response), ['patient_barcode', 'response_NR']]
patient_to_response = patient_to_response.drop_duplicates()
patient_to_response = patient_to_response.set_index('patient_barcode')['response_NR']

surv['response_NR'] = surv['patient_name'].replace(patient_to_response)
surv.loc[~surv['response_NR'].isin(['R', 'N']), 'response_NR'] = np.nan

surv['therapy_type'] = np.nan
surv.loc[surv['patient_name'].isin(ctla), 'therapy_type'] = 'anti-CTLA-4'

common_cols = surv.columns.intersection(data.columns)
data = data[common_cols]
surv = surv[common_cols]
print('Merging TIGER and TCGA data')
data['batch'] = 1
surv['batch'] = 2
data = pd.concat([data, surv])

# IMVigor210 Data
im = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/imvigor210/merged_log_fpkm.csv', index_col=0)
im['batch'] = 3

common_cols = im.columns.intersection(data.columns)
data = data[common_cols]
im = im[common_cols]

print('Merging TIGER, TCGA and IMVigor210 data')
data = pd.concat([data, im])
print(data.head())

#data.loc[~data['response_NR'].isin(['R', 'N']), 'response_NR'] = np.nan
data = data[data['response_NR'].isin(['R', 'N'])]

# Discard patients with many missing values
data = data[data.isna().sum(1) < 10000]

# Make sure patient names are unique
data.patient_name = data['dataset_id'].astype(str) + '-' + data.patient_name.astype(str)

clinical_cols = list(set(clinical_cols).intersection(common_cols)) 
clinical_cols.append('batch')
genes = [col for col in data.columns if col not in clinical_cols]

# Save data with also clinical cols and response
data.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/full_dataset_with_clinical.csv')

# Remove post-treatment samples
data = data[data['treatment_when'] != 'POST']

# Extract unique patients
patients = data[['patient_name', 'tissue', 'dataset_id']].drop_duplicates()

os.makedirs('/home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts', exist_ok=True)
os.makedirs('/home/lnemati/pathway_crosstalk/data/immunotherapy/tissues', exist_ok=True)

# Split patients based on tissue and dataset_id, cohorts have a minimum size
min_size = 20

# Get combinations of tissue and dataset_id that are actually present in the data
tissue_dataset_combinations = data[['tissue', 'dataset_id']].drop_duplicates()

def batch_correct(df, batch_col='batch_col', genes=genes):
    adata = AnnData(X = df[genes].fillna(np.log2(0.001)), obs = df[[batch_col]])
    corrected = combat(adata, key=batch_col, inplace=False)
    out = df.copy()
    out[genes] = corrected
    return out

#for tissue, dataset_id in tissue_dataset_combinations.itertuples(index=False):
#    print(tissue, dataset_id)
#    cohort_patients = patients.query('tissue == @tissue and dataset_id == @dataset_id')['patient_name']
#    # Check number of patients
#    if len(cohort_patients) < min_size:
#        print(f'Patients: {len(cohort_patients)}, skipping')
#        continue
#    # Check that both classes are present 
#    if not set(data.loc[data.patient_name.isin(cohort_patients), 'response_NR']).issuperset({'R', 'N'}):
#        print(f'Class not present, skipping')
#        continue
#   
#    cohort = data[data.patient_name.isin(cohort_patients)].copy()
#
#    # Only keep genes and response
#    cohort = cohort[genes + ['response_NR']]
#    print(tissue, dataset_id, f'n={len(cohort)}')
#    print('Nan values:', cohort.isna().sum().sum() / (len(cohort) * len(cohort.columns)))
#    print('Response:', cohort['response_NR'].value_counts())
#    print(cohort.head())
#
#    # Shuffle cohort
#    cohort = cohort.sample(frac=1, random_state=seed, ignore_index=False)
#
#    # Save cohort
#    filename = os.path.join('/home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts', f'{tissue}_{dataset_id}.csv'.lower().replace(' ', '_'))
#    cohort.to_csv(filename, index=True)
# Also save full dataset
#full_dataset = data[genes + ['response_NR']]
#full_dataset = full_dataset.sample(frac=1, random_state=seed)
#full_dataset.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts/full_dataset.csv')

# Also split based on tissue
for tissue in data['tissue'].unique():
    print(tissue)
    cohort_patients = patients.query('tissue == @tissue')['patient_name']
    # Check number of patients
    if len(cohort_patients) < min_size:
        print(f'Patients: {len(cohort_patients)}, skipping')
        continue
    # Check that both classes are present 
    if not set(data.loc[data.patient_name.isin(cohort_patients), 'response_NR']).issuperset({'R', 'N'}):
        print(f'Class not present')
        continue

    cohort = data[data.patient_name.isin(cohort_patients)].copy()
    
    # If more than one batch is present, perform batch correction
    if len(cohort['batch'].unique()) > 1 and BATCH_CORRECTION:
        print('Batches:', cohort['batch'].unique())
        cohort = batch_correct(cohort, batch_col='batch')

    # Only keep genes and response
    cohort = cohort[genes + ['response_NR']]
    print(tissue, f'n={len(cohort)}')
    print('Nan values:', cohort.isna().sum().sum() / (len(cohort) * len(cohort.columns)))
    print('Response:', cohort['response_NR'].value_counts())
    print(cohort.head())

    # Shuffle cohort
    cohort = cohort.sample(frac=1, random_state=seed, ignore_index=False)

    # Save cohort
    filename = os.path.join('/home/lnemati/pathway_crosstalk/data/immunotherapy/tissues', f'{tissue}.csv'.lower().replace(' ', '_'))
    cohort.to_csv(filename, index=True)

# Also save full dataset

# batch correction
if BATCH_CORRECTION:
    full_dataset = batch_correct(data, batch_col='batch')
full_dataset = data[genes + ['response_NR']]
full_dataset = full_dataset.sample(frac=1, random_state=seed)
full_dataset.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/cohorts/full_dataset.csv')
full_dataset.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/tissues/full_dataset.csv')

print('Done: split.py')
