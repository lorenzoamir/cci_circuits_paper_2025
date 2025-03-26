import os
import re
import pandas as pd
import numpy as np
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='Aggregate results')

parser.add_argument('--inputdir', type=str)

args = parser.parse_args()
parentdir = args.inputdir
outdir = os.path.join(parentdir, 'aggregated')
tissue = os.path.basename(parentdir)

print('Input directory:', parentdir)
print('Output directory:', outdir)
print('Tissue:', tissue)

os.makedirs(outdir, exist_ok=True)

dfs = []

possible_motifs = [
    '4_clique',
#    '4_no_crosstalk',
#    '4_triangle_extra',
#    '4_path',
#    '4_one_missing',
#    '4_cycle',
#    '3_path'
    '3_clique',
    'individual_ccis',
    'random_pairs'
]

for file in os.listdir(os.path.join(parentdir, 'metrics')):
    print(file)
    if not file.endswith('.csv'):
        continue
    if file.replace('.csv', '') not in possible_motifs:
        continue
    motif = file.replace('.csv', '')
    print(motif)
    df = pd.read_csv(os.path.join(parentdir, 'metrics', file), index_col=0)
    df['motif'] = motif
    print(df.auroc.mean())
    df['tissue'] = tissue
    dfs.append(df)
   
df = pd.concat(dfs)

# Get the genes in each motif
print('Getting genes')
df['all_genes'] = pd.Series(df.index).apply(lambda x: tuple(sorted(set(re.split('[&_+]', x))))).values

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
pairs['auprc1'] = ccc.loc[int1_genes, 'auprc'].values
pairs['auprc2'] = ccc.loc[int2_genes, 'auprc'].values

# Get difference between using both and max of only one
pairs['auroc_diff'] = pairs['auroc'] - pairs[['auroc1', 'auroc2']].max(axis=1)
pairs['auprc_diff'] = pairs['auprc'] - pairs[['auprc1', 'auprc2']].max(axis=1)

# Save the aggregated results
pairs = pairs.drop(columns=['all_genes'])
pairs = pairs.sort_values('auroc_diff', ascending=False)
# Reorder columns
cols = ['auroc', 'auroc1', 'auroc2', 'auprc', 'auprc1', 'auprc2', 'auroc_diff', 'auprc_diff', 'motif', 'tissue']
pairs = pairs[cols]
pairs.to_csv(os.path.join(outdir, 'aggregated.csv'))

print('Done: aggregate.py')
