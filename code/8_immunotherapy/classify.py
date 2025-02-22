import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
import re
import os
import argparse

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

train = train.drop(columns=['response_known'])

# response_NR and the clinical cols must be categorical and have the same categories in train and test
for col in clinical_cols:
    train[col] = train[col].astype('category')
    test[col] = test[col].astype('category').cat.set_categories(train[col].cat.categories)

# Find category number of 'R'
pos_label = list(train['response_NR'].cat.categories).index('R')

y_train = train['response_NR']
y_test  = test['response_NR']

# Read motif file
all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
# Get genes
all_motifs['all_genes'] = all_motifs.Interaction.apply(lambda x: re.split(r'[+_&]', x))

def scale_features(train, test, genes):
    # Always check that no patient is in both train and test
    assert set(train['patient_name']).isdisjoint(set(test['patient_name'])), "Error: Overlapping patients in train and test!"

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

def pca_dataset(train, test, genes, clinical_cols, n_pcs=0.9):
    # PCA
    pca = PCA(n_components=n_pcs, random_state=seed)
    X_train = pca.fit_transform(train[genes])
    X_test = pca.transform(test[genes])
    
    print('Using {} principal components'.format(pca.n_components_))

    # Convert back to DataFrame
    X_train = pd.DataFrame(X_train, index=train.index)
    X_test = pd.DataFrame(X_test, index=test.index)
    
    # TabPFN needs column names to be strings
    # Rename columns to PC1, PC2, ..., PCn
    X_train.columns = ['PC' + str(col) for col in X_train.columns]
    X_test.columns = ['PC' + str(col) for col in X_test.columns]

    return X_train, X_test

def train_test(train, test, genes, clinical_cols, n_pcs=0.9):
    X_train, X_test = pca_dataset(train, test, genes, clinical_cols, n_pcs)

    # Add back clinical columns
    X_train = X_train.join(train[clinical_cols])
    X_test = X_test.join(test[clinical_cols])

    categorical_feature_indices = list(range(X_train.shape[1]-len(clinical_cols), X_train.shape[1]))

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

if motif == 'whole_transcriptome':
    genes = set(train.columns)
    genes = list(genes - set(clinical_cols).union({'response_NR', 'patient_name', 'dataset_id'}))

    train, test = scale_features(train, test, genes)
    auroc, auprc = train_test(train, test, genes, clinical_cols, n_pcs=0.9)

    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')
    print()
    
    # Save results
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['whole_transcriptome'])
    results.to_csv('/home/lnemati/pathway_crosstalk/results/immunotherapy/whole_transcriptome.csv')

#if motif in all_motifs['Type'].unique():
#    motifdf = all_motifs[all_motifs['Type'] == motif]
#
#aurocs = {}
#
#for idx, row in motifdf.iterrows():
#    features = np.unique(row['Interaction'] + ['tumor_type'])
#
#    X_train = train[features]  
#    X_test  = test[features]
#
#    # Initialize a classifier
#    clf = TabPFNClassifier()
#    clf.fit(X_train, y_train)
#
#    # Predict probabilities
#    prediction_probabilities = clf.predict_proba(X_test)
#    auroc = roc_auc_score(y_test, prediction_probabilities[:, 1])
#
#    aurocs[idx] = auroc
#
#aurocs = pd.Series(aurocs).sort_values(ascending=False)
## make into a dataframe with two columns: interaction and auroc
#aurocs = pd.DataFrame(aurocs, columns=['auroc'])
#aurocs['interaction'] = aurocs.index
#aurocs = aurocs.reset_index(drop=True)
#aurocs = aurocs[['interaction', 'auroc']]
#aurocs['motif'] = motif
#
#outdir = '/home/lnemati/pathway_crosstalk/results/immunotherapy/aurocs'
#os.makedirs(outdir, exist_ok=True)
#
#aurocs.to_csv(os.path.join(outdir, f'{motif}_aurocs.csv'), index=False)
print('Done: classify.py')
