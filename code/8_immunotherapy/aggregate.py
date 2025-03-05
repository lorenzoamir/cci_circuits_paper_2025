import os
import re
import pandas as pd
import numpy as np
import argparse

parentdir = '/home/lnemati/pathway_crosstalk/results/immunotherapy/'

dfs = []
# Get aurocs
for file in os.listdir(os.path.join(parentdir, 'aurocs')):
    if not file.endswith('.csv'):
        continue
    if file == 'whole_transcriptome.csv':
        continue
    if file == 'all_ccis.csv':  
        continue
    motif = file.replace('.csv', '')
    print(motif)
    df = pd.read_csv(os.path.join(parentdir, 'aurocs', file), index_col=0)
    df['motif'] = motif
    dfs.append(df)
    
aurocs = pd.concat(dfs)
aurocs = aurocs[['mean', 'sem', 'motif']].rename(columns={'mean': 'auroc', 'sem': 'auroc_sem'})

dfs = []
for file in os.listdir(os.path.join(parentdir, 'auprcs')):
    if not file.endswith('.csv'):
        continue
    if file == 'whole_transcriptome.csv':
        continue
    if file == 'all_ccis.csv':
        continue
    motif = file.replace('.csv', '')
    print(motif)
    df = pd.read_csv(os.path.join(parentdir, 'auprcs', file), index_col=0)
    #df['motif'] = motif # it's already in the auroc df
    dfs.append(df)

auprcs = pd.concat(dfs)
auprcs = auprcs[['mean', 'sem']].rename(columns={'mean': 'auprc', 'sem': 'auprc_sem'})

# Join dataframes
print('Joining dataframes')
df = aurocs.join(auprcs)
print(df.head())

# Get the genes in each motif
print('Getting genes')
df['all_genes'] = pd.Series(df.index).apply(lambda x: tuple(sorted(set(re.split('[&_+]', x))))).values

print(df.motif.value_counts())

print(df.head())
ccc = df.query('motif == "individual_ccis"')
print(ccc.head())
pairs = df.query('motif != "individual_ccis"')
print(pairs.head())

print('Splitting interactions')
int1 = pd.Series(pairs.index).str.split('&', expand=True)[0]
int2 = pd.Series(pairs.index).str.split('&', expand=True)[1]

print('Inteactions 1:')
print(int1.head())
print('Inteactions 2:')
print(int2.head())

print('Getting genes in interactions')   
int1_genes = int1.apply(lambda x: tuple(sorted(set(re.split('[&_+]', x)))))
int2_genes = int2.apply(lambda x: tuple(sorted(set(re.split('[&_+]', x)))))
print(int1_genes.head())
print(int2_genes.head())

# Use all_genes as the index, sorted, no duplicates
ccc = ccc.set_index('all_genes')

print(ccc.head())
print(ccc.head())

pairs['auroc1'] = ccc.loc[int1_genes, 'auroc'].values
pairs['auroc2'] = ccc.loc[int2_genes, 'auroc'].values
pairs['auroc_sem1'] = ccc.loc[int1_genes, 'auroc_sem'].values
pairs['auroc_sem2'] = ccc.loc[int2_genes, 'auroc_sem'].values
pairs['auprc1'] = ccc.loc[int1_genes, 'auprc'].values
pairs['auprc2'] = ccc.loc[int2_genes, 'auprc'].values
pairs['auprc_sem1'] = ccc.loc[int1_genes, 'auprc_sem'].values
pairs['auprc_sem2'] = ccc.loc[int2_genes, 'auprc_sem'].values

# Get difference between using both and max of only one
pairs['auroc_diff'] = pairs['auroc'] - pairs[['auroc1', 'auroc2']].max(axis=1)
pairs['auroc_diff_err'] = np.sqrt(pairs['auroc_sem']**2 + pairs[['auroc_sem1', 'auroc_sem2']].max(axis=1)**2)
pairs['auprc_diff'] = pairs['auprc'] - pairs[['auprc1', 'auprc2']].max(axis=1)
pairs['auprc_diff_err'] = np.sqrt(pairs['auprc_sem']**2 + pairs[['auprc_sem1', 'auprc_sem2']].max(axis=1)**2)

# Save the aggregated results
outdir = os.path.join(parentdir, 'aggregated')
os.makedirs(outdir, exist_ok=True)
pairs = pairs.drop(columns=['all_genes'])
pairs = pairs.sort_values('auroc_diff', ascending=False)
# Reorder columns
cols = ['auroc', 'auroc1', 'auroc2', 'auprc', 'auprc1', 'auprc2', 'auroc_sem', 'auroc_sem1', 'auroc_sem2', 'auprc_sem', 'auprc_sem1', 'auprc_sem2', 'auroc_diff', 'auroc_diff_err', 'auprc_diff', 'auprc_diff_err', 'motif']
pairs = pairs[cols]
pairs.to_csv(os.path.join(outdir, 'aggregated.csv'))

print('Done: aggregate.py')
