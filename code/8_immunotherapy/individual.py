import pandas as pd
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.decomposition import PCA
from tabpfn import TabPFNClassifier
import signal
import random
import re
from ast import literal_eval
import os
import argparse

seed = 42

# Reproducibility
random.seed(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# Timeout setup
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Operation timed out.")

signal.signal(signal.SIGALRM, timeout_handler)
timeout_seconds = 20

n_folds = 10

parser = argparse.ArgumentParser()

parser.add_argument('--data', type=str, required=True)
parser.add_argument('--motif', type=str, required=True)
parser.add_argument('--outdir', type=str, required=False)

args = parser.parse_args()

results_dir = args.outdir
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

# Load data
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

# Get folds
cv = StratifiedKFold(n_splits=n_folds, shuffle=False)

if motif == 'individual_ccis':
    # Individual interactions 
    # Read cell-cell communication list
    ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
    ccc['all_genes'] = ccc['all_genes'].apply(literal_eval)

    aurocs = []
    auprcs = []
    probs = []

    ## DEBUG!
    #SUBSET = 2
    #SUBSET = min(SUBSET, ccc.shape[0])
    #print('Subsetting to ', SUBSET ,'random interactions')
    #ccc = ccc.sample(n=SUBSET, random_state=seed)

    for idx, row in ccc.iterrows():
        print('Processing interaction ', idx)

        genes = row['all_genes']
        
        if not set(genes).issubset(set(X.columns)):
            # skip if missing genes
            aurocs.append(np.nan)
            auprcs.append(np.nan)
            probs.append([np.nan] * X.shape[0])
            continue
        
        try:
            # Get prediction probabilities for each fold
            signal.alarm(timeout_seconds)
            prob = cross_val_predict(clf, X[genes], y, method='predict_proba', cv=cv)[:, 1]
            signal.alarm(0) # Cancel the alarm if successful

            auroc = roc_auc_score(y, prob, multi_class='ovr')
            auprc = average_precision_score(y, prob, pos_label=1)
        
        except TimeoutException:
            print('Timeout for interaction ', idx)
            prob = [np.nan] * X.shape[0] 
            auroc = np.nan
            auprc = np.nan

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

elif motif == 'random_pairs':
    random_pairs = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/random_pairs_of_interactions.csv', index_col=0)
    random_pairs['all_genes'] = random_pairs['all_genes'].apply(literal_eval)
    
    aurocs = []
    auprcs = []
    probs = []

    ## DEBUG!
    #SUBSET = 3
    #SUBSET = min(SUBSET, random_pairs.shape[0])
    #print('Subsetting to ', SUBSET ,'random pairs')
    #random_pairs = random_pairs.sample(n=SUBSET, random_state=seed)

    for idx, row in random_pairs.iterrows():
        genes = row['all_genes']

        if not set(genes).issubset(set(X.columns)):
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

    index = random_pairs.index

    # Save prediction probabilities (first row is target values)
    prediction_probabilities = pd.DataFrame(probs, index=index, columns=X.index)
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', 'random_pairs.csv'))

    # Save results
    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=index)
    results.to_csv(os.path.join(results_dir, 'metrics', 'random_pairs.csv'))

else:
    # Read motif file
    all_motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv')
    # Get genes
    all_motifs['all_genes'] = all_motifs.Interaction.apply(lambda x: re.split(r'[+_&]', x))

    if motif not in all_motifs['Type'].unique():
        raise ValueError(f'Motif {motif} not found in all_motifs')

    motifdf = all_motifs[all_motifs['Type'] == motif]

    # If there are more too many motifs, sample 1000
    if motifdf.shape[0] > 1000:
        print('Found ', motifdf.shape[0], ' motifs')
        print('Subsetting to 1000 random motifs')
        motifdf = motifdf.sample(n=1000, random_state=seed)

    aurocs = []
    auprcs = []
    probs = []

    ## DEBUG!
    #SUBSET = 2
    #SUBSET = min(SUBSET, motifdf.shape[0])
    #print('Subsetting to ', SUBSET ,'random motifs')
    #motifdf = motifdf.sample(n=SUBSET, random_state=seed)

    for idx, row in motifdf.iterrows():
        genes = row['all_genes']

        if not set(genes).issubset(set(X.columns)):
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

    # Save prediction probabilities (first row is target values)
    prediction_probabilities = pd.DataFrame(probs, index=motifdf['Interaction'].values, columns=X.index)
    prediction_probabilities = pd.concat([target, prediction_probabilities])
    prediction_probabilities.to_csv(os.path.join(results_dir, 'prediction_probabilities', f'{motif}.csv'))

    # Save results
    results = pd.DataFrame({'auroc': aurocs, 'auprc': auprcs}, index=motifdf['Interaction'].values) 
    results.to_csv(os.path.join(results_dir, 'metrics', f'{motif}.csv'))

print('Done: individual.py')
