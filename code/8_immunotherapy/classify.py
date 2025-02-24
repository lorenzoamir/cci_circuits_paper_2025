import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
from torch import set_float32_matmul_precision
import sys
import re
import os
import argparse

set_float32_matmul_precision('medium')

seed = 42
# Set seed for reproducibility
np.random.seed(seed)

parser = argparse.ArgumentParser()

parser.add_argument('--motif', type=str, required=True)

args = parser.parse_args()

motif = args.motif
print(motif)

clinical_cols = ['tissue', 'therapy_type', 'treatment_when']

dtypes = {
    'tissue': 'category',
    'therapy_type': 'category',
    'treatment_when': 'category',
    'response_NR': 'category',
}

# Load train and test set
train = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/train.csv', dtype=dtypes, index_col=0)
test = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/test.csv', dtype=dtypes, index_col=0) 

#train = train.drop(columns=['dataset_id', 'patient_name'])
#test = test.drop(columns=['dataset_id', 'patient_name'])

# response_NR and the clinical cols must be categorical and have the same categories in train and test
for col in clinical_cols:
    train[col] = train[col].astype('category')
    test[col] = test[col].astype('category').cat.set_categories(train[col].cat.categories)

# Find category number of 'R'
pos_label = list(train['response_NR'].cat.categories).index('R')

y_train = train['response_NR']
y_test  = test['response_NR']

train = train.drop(columns=['response_NR'])
test = test.drop(columns=['response_NR'])

def scale_features(train, test):
    # Always check that no patient is in both train and test

    genes = set(train.columns)
    genes = list(genes - set(clinical_cols).union({'response_NR', 'patient_name', 'dataset_id'}))

    # Scaling
    scaler = StandardScaler()
    X_train = scaler.fit_transform(train[genes])
    X_test = scaler.transform(test[genes])

    # Convert back to DataFrame
    X_train = pd.DataFrame(X_train, index=train.index, columns=genes)
    X_test = pd.DataFrame(X_test, index=test.index, columns=genes)
  
    # Add back clinical columns
    train = X_train.join(train[clinical_cols])
    test = X_test.join(test[clinical_cols])

    return train, test

def pca_dataset(train, test, genes, clinical_cols, n_pcs=30):
    # PCA
    pca = PCA(n_components=n_pcs, random_state=seed)
    X_train = pca.fit_transform(train[genes])
    X_test = pca.transform(test[genes])
    
    # Convert back to DataFrame
    X_train = pd.DataFrame(X_train, index=train.index)
    X_test = pd.DataFrame(X_test, index=test.index)
    
    # TabPFN needs column names to be strings
    # Rename columns to PC1, PC2, ..., PCn
    X_train.columns = ['PC' + str(col) for col in X_train.columns]
    X_test.columns = ['PC' + str(col) for col in X_test.columns]

    return X_train, X_test

def train_test(train, test, genes, clinical_cols, n_pcs=None):
    if n_pcs is not None:
        X_train, X_test = pca_dataset(train, test, genes, clinical_cols, n_pcs)
        # Add back clinical columns
        X_train = X_train.join(train[clinical_cols])
        X_test = X_test.join(test[clinical_cols])
    else:
        X_train = train[genes + clinical_cols]
        X_test = test[genes + clinical_cols]

    categorical_feature_indices = list(range(X_train.shape[1]-len(clinical_cols), X_train.shape[1]))

    # DEBUG:
    print('X_train:')
    print(X_train.head())
    print('X_test:')
    print(X_test.head())
    print('columns:', X_train.columns)
    print('categorical_feature_indices:', categorical_feature_indices)
    print('categorical_columns:', X_train.columns[categorical_feature_indices])

    clf = TabPFNClassifier(
        random_state=seed,
        balance_probabilities='True',
        categorical_features_indices=categorical_feature_indices,
        device='cuda',
    )

    clf.fit(X_train, y_train)

    # Predict probabilities
    prediction_probabilities = clf.predict_proba(X_test)

    # Compute AUROC
    auroc = roc_auc_score(y_test, prediction_probabilities[:, 1])

    # Compute AUPRC
    auprc = average_precision_score(y_test, prediction_probabilities[:, 1], pos_label='R') 

    return auroc, auprc

train, test = scale_features(train, test)

if motif == 'whole_transcriptome':
    genes = set(train.columns)
    genes = list(genes - set(clinical_cols).union({'response_NR', 'patient_name', 'dataset_id'}))

    # Perform data augmentation by randomly setting genes to 0 with a probability of 0.1
    fake_train = train.copy()
    fake_y_train = y_train.copy()
    mask = np.random.choice([0, 1], size=(fake_train.shape[0], len(genes)), p=[0.1, 0.9])

    # Set cells to 0
    fake_train[genes] = fake_train[genes].values * mask
    
    # Train on both real and fake data
    train = pd.concat([train, fake_train])
    y_train = pd.Series(np.concatenate([y_train.values, fake_y_train.values]), index=train.index)

    auroc, auprc = train_test(train, test, genes, clinical_cols, n_pcs=30)

    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')
    print()
    
    # Save results
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['whole_transcriptome'])
    results.to_csv('/home/lnemati/pathway_crosstalk/results/immunotherapy/whole_transcriptome.csv')

else:
    # Read motif file
    all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
    # Get genes
    all_motifs['all_genes'] = all_motifs.Interaction.apply(lambda x: re.split(r'[+_&]', x))

    if motif not in all_motifs['Type'].unique():
        raise ValueError(f'Motif {motif} not found in all_motifs')

    motifdf = all_motifs[all_motifs['Type'] == motif]

    aurocs = []
    auprcs = []

    # DEBUG!
    SUBSET = 1000
    SUBSET = min(SUBSET, motifdf.shape[0])
    print('Subsetting to ', SUBSET ,'random motifs', file=sys.stdout)
    print('Subsetting to ', SUBSET ,'random motifs', file=sys.stderr)

    motifdf = motifdf.sample(n=SUBSET, random_state=seed)

    for idx, row in motifdf.iterrows():
        genes = row['all_genes']

        auroc, auprc = train_test(train, test, genes, clinical_cols, n_pcs=None)
        
        aurocs.append(auroc)
        auprcs.append(auprc)

    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=motifdf['Interaction'].values)
    results.to_csv(f'/home/lnemati/pathway_crosstalk/results/immunotherapy/{motif}.csv')

print('Done: classify.py')
