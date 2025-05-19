import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.decomposition import PCA
from tabpfn import TabPFNClassifier
import random
import re
import pickle
from ast import literal_eval
import os
import argparse

seed = 42
# Set seed for reproducibility
random.seed(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# n_folds = 15 # Default number of folds
n_pcs = 0.95

parser = argparse.ArgumentParser()

parser.add_argument('--data', type=str, required=True)
parser.add_argument('--features', type=str, required=True)
parser.add_argument('--outdir', type=str, required=False)

args = parser.parse_args()

results_dir = args.outdir
cohort_name = args.data.split('/')[-1].split('.')[0]
print(cohort_name)
results_dir = os.path.join(results_dir, cohort_name)
os.makedirs(results_dir, exist_ok=True)
os.makedirs(os.path.join(results_dir, 'prediction_probabilities'), exist_ok=True)
os.makedirs(os.path.join(results_dir, 'metrics'), exist_ok=True)
os.makedirs(os.path.join(results_dir, 'models'), exist_ok=True)

features = args.features
print(features)

dtypes = {
    'response_NR': 'category',
}

# Load train and test set
X = pd.read_csv(args.data, index_col=0, dtype=dtypes)
target = X['response_NR']
y = np.where(target == 'R', 1, 0)
X = X.drop(columns=['response_NR'])

# Drop columns where all values are nans or the same
X = X.dropna(axis=1, how='all')
X = X.loc[:, X.nunique() != 1]

# Make target a dataframe of a single row, with X.index as columns
target.name = 'target'
target = target.to_frame().T

# Find category number of 'R'
n_samples = y.shape[0]

# Init classifier
clf = TabPFNClassifier(
    random_state=seed,
    balance_probabilities='True',
    device='cuda',
)

# One fold for each sample in the minority class
n_folds = min(np.sum(y == 1), np.sum(y == 0)) 

# Get folds
cv = StratifiedKFold(n_splits=n_folds, shuffle=False)

def pca_dataset(X, n_pcs=None):
    # Default is to keep all PCs, if n_pcs is a fraction, keep PCs that explain that fraction of the variance
    # PCA
    pca = PCA(random_state=seed)
    X_pca = pca.fit_transform(X.fillna(np.log2(0.001)))
    total_pcs = pca.n_components_ 

    # Convert back to DataFrame
    X_pca = pd.DataFrame(X_pca, index=X.index)
     
    # If n_pcs is a fraction, keep PCs that explain that fraction of the variance
    if n_pcs is None:
        n_pcs = total_pcs
    elif n_pcs < 1:
        explained_variance_ratio_cumsum = np.cumsum(pca.explained_variance_ratio_)
        n_pcs = np.argmax(explained_variance_ratio_cumsum > n_pcs) + 1

    # TabPFN has a max of 500 features
    n_pcs = min(n_pcs, 500)
    X_pca = X_pca.iloc[:, :n_pcs]
            
    print('Number of PCs:', pca.n_components_)
    print('Explained variance:', np.sum(pca.explained_variance_ratio_))

    # TabPFN needs column names to be strings
    # Rename columns to PC1, PC2, ..., PCn
    X_pca.columns = ['PC' + str(col) for col in X_pca.columns]
    
    return X_pca

if features == 'whole_transcriptome':
    X = pca_dataset(X, n_pcs=n_pcs)
    
    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]
    auroc = roc_auc_score(y, probs, multi_class='ovr')
    auprc = average_precision_score(y, probs, pos_label=1)
    
    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')

    # Save prediction probabilities (first row is target values)
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['whole_transcriptome']).T
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'whole_transcriptome.csv')) 
    
    # Save metrics
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['whole_transcriptome'])
    results.to_csv(os.path.join(results_dir, 'metrics', 'whole_transcriptome.csv'))

    # Save model
    with open(os.path.join(results_dir, 'models', 'whole_transcriptome.pkl'), 'wb') as f:
        pickle.dump(clf, f)

elif features == 'all_ccis':
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    genes = set(ccc['all_genes'].sum())
    genes = list(genes.intersection(set(X.columns)))
    
    X = pca_dataset(X[genes], n_pcs=n_pcs)

    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]
    auroc = roc_auc_score(y, probs, multi_class='ovr')
    auprc = average_precision_score(y, probs, pos_label=1)
    
    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')

    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['all_ccis']).T
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'all_ccis.csv'))

    # Save metrics
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['all_ccis'])
    results.to_csv(os.path.join(results_dir, 'metrics', 'all_ccis.csv'))

    # Save model
    with open(os.path.join(results_dir, 'models', 'all_ccis.pkl'), 'wb') as f:
        pickle.dump(clf, f)

elif features == 'all_cliques':
    all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
    cliques = all_motifs[all_motifs['Type'].isin(['3_clique', '4_clique'])]
    cliques['all_genes'] = cliques.Interaction.apply(lambda x: re.split(r'[+_&]', x))
    genes = set(cliques['all_genes'].sum())
    genes = list(genes.intersection(set(X.columns)))

    #X = pca_dataset(X[genes], n_pcs=n_pcs)
    X = X[genes]

    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]
    auroc = roc_auc_score(y, probs, multi_class='ovr')
    auprc = average_precision_score(y, probs, pos_label=1)

    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')

    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['all_cliques']).T
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'all_cliques.csv'))

    # Save metrics
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['all_cliques'])
    results.to_csv(os.path.join(results_dir, 'metrics', 'all_cliques.csv'))
    
    # Save model
    with open(os.path.join(results_dir, 'models', 'all_cliques.pkl'), 'wb') as f:
        pickle.dump(clf, f)

elif features == 'signatures':
    # Use genes from TIGER signatures
    signatures = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/tiger_signatures/all_signatures.csv', index_col=0)
    signatures['all_genes'] = signatures['genes'].apply(literal_eval)

    for signature_idx in signatures.index:
        signature = signature_idx.lower().replace(' ', '_')

        genes = set(signatures.loc[signature_idx, 'all_genes'])
        genes = list(genes.intersection(set(X.columns)))

        X_signature = X[genes]

        # If there are more than 500 genes, use PCA
        if len(genes) > 500:
            X_signature = pca_dataset(X_signature, n_pcs=n_pcs)
            signature = f'{signature}_pca'
        
        probs = cross_val_predict(clf, X_signature, y, method='predict_proba', cv=cv)[:, 1]
        auroc = roc_auc_score(y, probs, multi_class='ovr')
        auprc = average_precision_score(y, probs, pos_label=1)

        print(f'{signature} AUROC: {auroc}')
        print(f'{signature} AUPRC: {auprc}')

        # Save prediction probabilities
        prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=[signature]).T
        prediction_probabilities = pd.concat([target, prediction_probabilities])
        prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', f'{signature}.csv'))

        # Save metrics
        results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=[signature])
        results.to_csv(os.path.join(results_dir, 'metrics', f'{signature}.csv'))
        
        # Save model
        with open(os.path.join(results_dir, 'models', f'{signature}.pkl'), 'wb') as f:
            pickle.dump(clf, f)

print('Done: features.py')
