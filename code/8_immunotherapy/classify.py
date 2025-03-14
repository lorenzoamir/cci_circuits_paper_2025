import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.decomposition import PCA
from tabpfn import TabPFNClassifier
import re
from ast import literal_eval
import os
import argparse

seed = 42
# Set seed for reproducibility
np.random.seed(seed)

n_folds = 10

parser = argparse.ArgumentParser()

parser.add_argument('--data', type=str, required=True)
parser.add_argument('--motif', type=str, required=True)

args = parser.parse_args()

results_dir = '/home/lnemati/pathway_crosstalk/results/immunotherapy/cohorts/'
cohort_name = args.data.split('/')[-1].split('.')[0]
print(cohort_name)
results_dir = os.path.join(results_dir, cohort_name)
os.makedirs(results_dir, exist_ok=True)
os.makedirs(os.path.join(results_dir, 'prediction_probabilities'), exist_ok=True)
os.makedirs(os.path.join(results_dir, 'metrics'), exist_ok=True)

motif = args.motif
print(motif)

dtypes = {
    'response_NR': 'category',
}

# Load train and test set
X = pd.read_csv(args.data, index_col=0, dtype=dtypes)
target = X['response_NR']
y = np.where(target == 'R', 1, 0)
X = X.drop(columns=['response_NR'])

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

# Get folds
cv = StratifiedKFold(n_splits=n_folds, shuffle=False)

def pca_dataset(X, n_pcs=None):
    # Default number of PCs is the number of samples
    # but if those explain more than 90% of the variance, use that number
    if n_pcs is None:
        n_pcs = min(X.shape[0], X.shape[1])

    # PCA
    print('Trying wih n_pcs:', n_pcs)
    pca = PCA(n_components=n_pcs, random_state=seed)
    X_pca = pca.fit_transform(X.fillna(np.log2(0.001)))

    # If n_pcs explain more than 90% of the variance, only keep those that explain 90%
    if n_pcs > 1:
        explained_variance_ratio = pca.explained_variance_ratio_
        explained_variance_ratio_cumsum = np.cumsum(explained_variance_ratio)
        total_variance = explained_variance_ratio_cumsum[-1]
        if total_variance > 0.9: 
            # Set n_pcs to the number of PCs that explain 90% of the variance
            n_pcs = np.argmax(explained_variance_ratio_cumsum > 0.9) + 1
            pca = PCA(n_components=n_pcs, random_state=seed)
            X_pca = pca.fit_transform(X.fillna(np.log2(0.001)))
            
    print('Number of PCs:', pca.n_components_)
    print('Explained variance:', np.sum(pca.explained_variance_ratio_))

    # Convert back to DataFrame
    X_pca = pd.DataFrame(X_pca, index=X.index)
    
    # TabPFN needs column names to be strings
    # Rename columns to PC1, PC2, ..., PCn
    X_pca.columns = ['PC' + str(col) for col in X_pca.columns]
    
    return X_pca

if motif == 'whole_transcriptome':
    X = pca_dataset(X)
    
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

elif motif == 'all_ccis':
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    genes = set(ccc['all_genes'].sum())
    genes = list(genes.intersection(set(X.columns)))
    
    X = pca_dataset(X[genes])

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

elif motif == 'all_motifs':
    # Read motifs file
    all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
    # Get all genes that appear in any motif
    all_motifs['all_genes'] = all_motifs.Interaction.apply(lambda x: re.split(r'[+_&]', x))
    genes = set(all_motifs['all_genes'].sum())
    genes = list(genes.intersection(set(X.columns)))

    X = pca_dataset(X[genes])

    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]
    auroc = roc_auc_score(y, probs, multi_class='ovr')
    auprc = average_precision_score(y, probs, pos_label=1)

    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')

    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['all_motifs']).T
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'all_motifs.csv'))

    # Save metrics
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['all_motifs'])
    results.to_csv(os.path.join(results_dir, 'metrics', 'all_motifs.csv'))

elif motif == 'individual_ccis':
    # Individual interactions 
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    ## DEBUG!
    #SUBSET = 20
    #SUBSET = min(SUBSET, ccc.shape[0])
    #print('Subsetting to ', SUBSET ,'random interactions')
    #ccc = ccc.sample(n=SUBSET, random_state=seed)

    aurocs = []
    auprcs = []
    probs = []

    for idx, row in ccc.iterrows():
        genes = row['all_genes']
        
        if not set(genes).issubset(set(X.columns)):
            # skip if missing genes
            aurocs.append(np.nan)
            auprcs.append(np.nan)
            probs.append([np.nan] * X.shape[0])
            continue

        # Get prediction probabilities for each fold
        prob = cross_val_predict(clf, X[genes], y, method='predict_proba', cv=cv)[:, 1]
        auroc = roc_auc_score(y, prob, multi_class='ovr')
        auprc = average_precision_score(y, prob, pos_label=1)

        aurocs.append(auroc)
        auprcs.append(auprc)
        probs.append(prob)
    
    index = ccc['interaction'].str.replace('_', 'TEMP_REPLACE').str.replace('+', '_').str.replace('TEMP_REPLACE', '+').values
    
    # Save prediction probabilities (first row is target values)
    prediction_probabilities = pd.DataFrame(probs, index=index, columns=X.index)
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'individual_ccis.csv'))

    #Save results
    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=index)
    results.to_csv(os.path.join(results_dir, 'metrics', 'individual_ccis.csv'))
    
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
    probs = []

    # DEBUG!
    #SUBSET = 1000
    #SUBSET = min(SUBSET, motifdf.shape[0])
    #print('Subsetting to ', SUBSET ,'random motifs')

    #motifdf = motifdf.sample(n=SUBSET, random_state=seed)

    for idx, row in motifdf.iterrows():
        genes = row['all_genes']

        if not set(genes).issubset(set(train.columns)):
            aurocs.append(np.nan)
            auprcs.append(np.nan)
            probs.append([np.nan] * test.shape[0])
            continue
        
        X_train = train[genes + clinical_cols]
        X_test = test[genes + clinical_cols]
        auroc, auprc, prediction_probabilities, clf = train_test(X_train, X_test, genes, clinical_cols)
       
        aurocs.append(auroc)
        auprcs.append(auprc)
        probs.append(prediction_probabilities[:, 1])

    # Save results
    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=motifdf['Interaction'].values)
    results.to_csv(os.path.join(results_dir, f'{motif}.csv'))
    
    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=motifdf['Interaction'].values, columns=test.index)
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', f'{motif}.csv'))

print('Done: classify.py')
