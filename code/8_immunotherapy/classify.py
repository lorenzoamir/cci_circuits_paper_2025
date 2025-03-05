import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from tabpfn import TabPFNClassifier
from torch import set_float32_matmul_precision, manual_seed
from ast import literal_eval
from scipy.stats import sem
import sys
import re
import os
import random
import argparse

#set_float32_matmul_precision('medium')

# Set seed for reproducibility
seed = 42
random.seed(seed)
np.random.seed(seed)
manual_seed(seed)

parser = argparse.ArgumentParser()

parser.add_argument('--motif', type=str, required=True)
args = parser.parse_args()

results_dir = '/home/lnemati/pathway_crosstalk/results/immunotherapy'
os.makedirs(results_dir, exist_ok=True)
os.makedirs(os.path.join(results_dir, 'prediction_probabilities'), exist_ok=True)
os.makedirs(os.path.join(results_dir, 'aurocs'), exist_ok=True)
os.makedirs(os.path.join(results_dir, 'auprcs'), exist_ok=True)

motif = args.motif
print(motif)

#clinical_cols = ['tissue', 'therapy_type', 'treatment_when']
categorical_cols = ['patient_name', 'dataset_id', 'response_NR', 'tissue', 'therapy_type', 'treatment_when']

clinical_cols = ['tissue', 'therapy_type', 'treatment_when']
cols_to_remove = list(set(categorical_cols) - set(clinical_cols))

dtypes = {col: 'category' for col in categorical_cols}

# Load train and test set
#data = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/batch_corrected.csv', dtype=dtypes, index_col=0)
data = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/batch_corrected.csv', index_col=0, dtype=dtypes)
splits = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/splits.csv', index_col=0)

# Make sure data index is made of strings
data.index = data.index.astype(str)

# response_NR and the clinical cols must be categorical and have the same categories in train and test
for col in clinical_cols:
    data[col] = data[col].astype('category')

# Find category number of 'R'
pos_label = list(data['response_NR'].cat.categories).index('R')
#n_pcs = 30
n_pcs = 0.9

def get_train_test(data, splits, split_column):
    # splits has the same index as data and columns 'split0', 'split1', ...,
    # with values 'train' or 'test', assign samples to train or test set
  
    train_idx = list(splits.index[splits[split_column] == 'train'])
    test_idx = list(splits.index[splits[split_column] == 'test'])
    
    train = data.loc[train_idx]
    test = data.loc[test_idx]

    y_train = train['response_NR']
    y_test = test['response_NR']

    train = train.drop(columns=['response_NR'])
    test = test.drop(columns=['response_NR'])

    return train, test, y_train, y_test

def scale_features(train, test):
    # Always check that no patient is in both train and test

    genes = set(train.columns)
    genes = list(genes - set(clinical_cols).union(cols_to_remove))

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
    # Fill nans with 0
    train[genes] = train[genes].fillna(0)
    test[genes] = test[genes].fillna(0)

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

    return auroc, auprc, prediction_probabilities

if motif == 'whole_transcriptome':
    aurocs = pd.DataFrame(index=['whole_transcriptome'], columns=splits.columns)
    auprcs = pd.DataFrame(index=['whole_transcriptome'], columns=splits.columns)
    
    genes = set(data.columns)
    genes = list(genes - set(clinical_cols).union(cols_to_remove))

    for splitn in splits.columns:
        train, test, y_train, y_test = get_train_test(data, splits, splitn)
        train, test = scale_features(train, test)

        auroc, auprc, prediction_probabilities = train_test(train, test, genes, clinical_cols, n_pcs=n_pcs)
        
        # Add auroc and auprc to the dataframes
        aurocs.at['whole_transcriptome', splitn] = auroc
        auprcs.at['whole_transcriptome', splitn] = auprc

elif motif == 'all_ccis':
    aurocs = pd.DataFrame(index=['all_ccis'], columns=splits.columns) 
    auprcs = pd.DataFrame(index=['all_ccis'], columns=splits.columns)

    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)
    
    genes = set(ccc['all_genes'].sum())
    genes = list(genes.intersection(set(data.columns)))

    for splitn in splits.columns:
        train, test, y_train, y_test = get_train_test(data, splits, splitn)
        train, test = scale_features(train, test)

        auroc, auprc, prediction_probabilities = train_test(train, test, genes, clinical_cols, n_pcs=n_pcs)

        # Add auroc and auprc to the dataframes
        aurocs.at['all_ccis', splitn] = auroc
        auprcs.at['all_ccis', splitn] = auprc

elif motif == 'individual_ccis':
    # Individual interactions 
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)

    ## DEBUG
    #SUBSET = 20
    #SUBSET = min(SUBSET, ccc.shape[0])
    #print('Subsetting to ', SUBSET ,'random interactions', file=sys.stdout)
    #print('Subsetting to ', SUBSET ,'random interactions', file=sys.stderr)
    #ccc = ccc.sample(n=SUBSET, random_state=seed)

    aurocs = pd.DataFrame(index=ccc['interaction'], columns=splits.columns)
    auprcs = pd.DataFrame(index=ccc['interaction'], columns=splits.columns)

    for splitn in splits.columns:
        train, test, y_train, y_test = get_train_test(data, splits, splitn)
        train, test = scale_features(train, test)

        for idx, row in ccc.iterrows():
            genes = row['all_genes']
            
            if not set(genes).issubset(set(train.columns)):
                aurocs.at[row['interaction'], splitn] = np.nan
                auprcs.at[row['interaction'], splitn] = np.nan
                continue

            auroc, auprc, prediction_probabilities = train_test(train, test, genes, clinical_cols, n_pcs=None)

            aurocs.at[row['interaction'], splitn] = auroc
            auprcs.at[row['interaction'], splitn] = auprc

    # Replace '+' and '_' in the index
    index = ccc['interaction'].str.replace('_', 'TEMP_REPLACE').str.replace('+', '_').str.replace('TEMP_REPLACE', '+').values
    aurocs.index = index
    auprcs.index = index

else:
    # Read motif file
    all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
    # Get genes
    all_motifs['all_genes'] = all_motifs.Interaction.apply(lambda x: re.split(r'[+_&]', x))

    if motif not in all_motifs['Type'].unique():
        raise ValueError(f'Motif {motif} not found in all_motifs')

    motifdf = all_motifs[all_motifs['Type'] == motif]

    ## DEBUG!
    #SUBSET = 20
    #SUBSET = min(SUBSET, motifdf.shape[0])
    #print('Subsetting to ', SUBSET ,'random motifs', file=sys.stdout)
    #print('Subsetting to ', SUBSET ,'random motifs', file=sys.stderr)
    #motifdf = motifdf.sample(n=SUBSET, random_state=seed)

    aurocs = pd.DataFrame(index=motifdf['Interaction'].values, columns=splits.columns)
    auprcs = pd.DataFrame(index=motifdf['Interaction'].values, columns=splits.columns)
    
    for splitn in splits.columns:
        train, test, y_train, y_test = get_train_test(data, splits, splitn)
        train, test = scale_features(train, test)

        for idx, row in motifdf.iterrows():
            genes = row['all_genes']

            auroc, auprc, prediction_probabilities = train_test(train, test, genes, clinical_cols, n_pcs=None)
           
            aurocs.at[row['Interaction'], splitn] = auroc
            auprcs.at[row['Interaction'], splitn] = auprc

# Get mean, std and standard error of the mean
splits_cols = [col for col in splits.columns if col.startswith('split')]
aurocs['mean'] = aurocs[splits_cols].mean(axis=1)
aurocs['std'] = aurocs[splits_cols].std(axis=1)
aurocs['sem'] = aurocs[splits_cols].apply(sem, axis=1)

auprcs['mean'] = auprcs[splits_cols].mean(axis=1)
auprcs['std'] = auprcs[splits_cols].std(axis=1)
auprcs['sem'] = auprcs[splits_cols].apply(sem, axis=1)

# Sort by mean auroc, reorder columns
aurocs = aurocs.sort_values(by='mean', ascending=False)[['mean', 'std', 'sem'] + splits_cols] # sort
auprcs = auprcs.loc[aurocs.index, ['mean', 'std', 'sem'] + splits_cols] # match order of aurocs

# Save results
aurocs.to_csv(os.path.join(results_dir, 'aurocs', f'{motif}.csv'))
auprcs.to_csv(os.path.join(results_dir, 'auprcs', f'{motif}.csv'))

print('Done: classify.py')
