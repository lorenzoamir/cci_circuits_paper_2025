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

motif = args.motif
print(motif)

dtypes = {
    'response_NR': 'category',
}

# Load train and test set
X = pd.read_csv(args.data, index_col=0, dtype=dtypes)
y = X['response_NR']
X = X.drop(columns=['response_NR'])

# Convert 'N' and 'R' to 0 and 1
y = np.where(y == 'R', 1, 0)

# Find category number of 'R'
n_pcs = 0.9

# Init classifier
clf = TabPFNClassifier(
    random_state=seed,
    balance_probabilities='True',
    device='cuda',
)

# Get folds
cv = StratifiedKFold(n_splits=n_folds, shuffle=False)

def pca_dataset(X, n_pcs=0.9):
    # PCA
    pca = PCA(n_components=n_pcs, random_state=seed)
    X_pca = pca.fit_transform(X.fillna(np.log2(0.001)))
    
    print('Number of PCs:', pca.n_components_)

    # Convert back to DataFrame
    X_pca = pd.DataFrame(X_pca, index=X.index)
    
    # TabPFN needs column names to be strings
    # Rename columns to PC1, PC2, ..., PCn
    X_pca.columns = ['PC' + str(col) for col in X_pca.columns]
    
    return X_pca

if motif == 'whole_transcriptome':
    X = pca_dataset(X, n_pcs=n_pcs)
    
    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]
    # Save prediction probabilities 
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['whole_transcriptome']).T
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'whole_transcriptome.csv')) 
    
    # Compute AUROC and AUPRC
    auroc = roc_auc_score(y, probs)
    auprc = average_precision_score(y, probs, pos_label=1)
    
    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')
    print()
    
    # Save results
    results = pd.DataFrame({'auroc': [auroc], 'auprc': [auprc]}, index=['whole_transcriptome'])
    
    # Print real values vs prediction probabilities
    df = pd.DataFrame({'real': y, 'prediction': probs}, index=X.index)
    print(df)

elif motif == 'all_ccis':
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    genes = set(ccc['all_genes'].sum())
    genes = list(genes.intersection(set(X.columns)))
    
    X = pca_dataset(X[genes], n_pcs=n_pcs)

    # Get prediction probabilities for each fold
    probs = cross_val_predict(clf, X, y, method='predict_proba', cv=cv)[:, 1]

    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=X.index, columns=['all_ccis']).T
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'all_ccis.csv'))

    # Compute AUROC and AUPRC
    auroc = roc_auc_score(y, probs)
    auprc = average_precision_score(y, probs, pos_label=1)
    
    print(f'AUROC: {auroc}')
    print(f'AUPRC: {auprc}')
    print()

elif motif == 'individual_ccis':
    # Individual interactions 
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    aurocs = []
    auprcs = []
    probs = []

    for idx, row in ccc.iterrows():
        genes = row['all_genes']
        
        if not set(genes).issubset(set(train.columns)):
            aurocs.append(np.nan)
            auprcs.append(np.nan)
            probs.append([np.nan] * test.shape[0])
            continue
        
        X_train = train[genes + clinical_cols]
        X_test = test[genes + clinical_cols]
        auroc, auprc, prediction_probabilities, clf = train_test(X_train, X_test, genes, clinical_cols, n_pcs=None)

        aurocs.append(auroc)
        auprcs.append(auprc)
        probs.append(prediction_probabilities[:, 1])
    
    index = ccc['interaction'].str.replace('_', 'TEMP_REPLACE').str.replace('+', '_').str.replace('TEMP_REPLACE', '+').values
    
    #Save results
    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=index)
    results.to_csv(os.path.join(results_dir, 'individual_ccis.csv'))
    
    # Save prediction probabilities
    prediction_probabilities = pd.DataFrame(probs, index=index, columns=test.index)
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'individual_ccis.csv'))

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
    #print('Subsetting to ', SUBSET ,'random motifs', file=sys.stdout)
    #print('Subsetting to ', SUBSET ,'random motifs', file=sys.stderr)

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
        auroc, auprc, prediction_probabilities, clf = train_test(X_train, X_test, genes, clinical_cols, n_pcs=None)
       
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
