import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare hub genes from WGCNAs')

parser.add_argument('--dir_list', type=str, help='List of directories to compare, delimited by commas')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results'

perc = 0.05

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]

check_files = [
    'degree_df.csv.gz',
    'interactions.csv',
]

# Sanity check
for d in dir_list:
    if not os.path.exists(d):
        # error message to standard error and standard output
        print('Error: Directory {} does not exist'.format(d), file=sys.stderr)
        print('Error: Directory {} does not exist'.format(d), file=sys.stdout)
    # Check if all files are present
    for f in check_files:
        if not os.path.exists(os.path.join(d, f)):
            print('Error: File {} does not exist in directory {}'.format(f, d), file=sys.stderr)
            print('Error: File {} does not exist in directory {}'.format(f, d), file=sys.stdout)

def get_name_condition(d):
    if '/tumor/' in d: 
        condition = 'tumor'
    elif '/normal/' in d:
        condition = 'normal'
    else:
        print('Error: Directory {} does not contain /tumor or /normal'.format(d), file=sys.stderr)
        print('Error: Directory {} does not contain /tumor or /normal'.format(d), file=sys.stdout)
        return None, None
    
    name = d.split('/')[-1]
    return name, condition

global_hubs = {}
all_genes = set()
all_interactors = set()

normal_tissues = []
tumor_tissues = []

print('Reading all files')
for d in dir_list:

    print(d)
    name, condition = get_name_condition(d)
    print(name)
   
    if condition == 'normal':
        normal_tissues.append(name)
    elif condition == 'tumor':
        tumor_tissues.append(name)
    
    degree_df = pd.read_csv(os.path.join(d, 'degree_df.csv.gz'), index_col=0)
    interactions = pd.read_csv(os.path.join(d, 'interactions.csv'))

    # Add all genes to set
    all_genes.update(set(degree_df.index))
    
    # Create set of all interactors by taking all unique values that appear in 
    # columns interactor1, interactor2, ..., interactor7
    for i in range(1, 8):
        all_interactors.update(set(interactions['interactor{}'.format(i)]))

    print('Number of interactors: {}'.format(len(all_interactors)))

    # Get top 5% of genes by degree
    print('Getting top 5% of genes by degree')
    global_hubs[name] = degree_df.sort_values('degree', ascending=False).head(int(perc * len(degree_df))).index
    print(global_hubs[name])

print('Done reading all files')

print('Scoring each gene')
# Score each gene by how many times it appears in the top 5% of genes
genes_scores = pd.DataFrame(index=list(all_genes), columns=['global_hub_tumor', 'global_hub_normal', 'interactor'])

for gene in all_genes:
    global_hub_tumor = 0
    global_hub_normal = 0
    genes_scores.loc[gene, 'interactor'] = 1 if gene in all_interactors else 0

    for name in global_hubs:
        if gene in global_hubs[name]:
            if name in tumor_tissues:
                global_hub_tumor += 1
            elif name in normal_tissues:
                global_hub_normal += 1

    genes_scores.loc[gene, 'global_hub_tumor'] = global_hub_tumor
    genes_scores.loc[gene, 'global_hub_normal'] = global_hub_normal

# Use freq instad of count
genes_scores['global_hub_tumor'] = genes_scores['global_hub_tumor'] / len(tumor_tissues)
genes_scores['global_hub_normal'] = genes_scores['global_hub_normal'] / len(normal_tissues)

print('Done scoring each gene')

print('Removing genes that are never hubs')

# Remove genes that are never hubs
genes_scores = genes_scores[(genes_scores['global_hub_tumor'] > 0) | (genes_scores['global_hub_normal'] > 0)]
   
# Sort by differece between tumor and normal
genes_scores['diff'] = genes_scores['global_hub_tumor'] - genes_scores['global_hub_normal']
genes_scores = genes_scores.sort_values('diff', ascending=False)
genes_scores.drop(columns='diff', inplace=True)

# Save results
genes_scores.to_csv(os.path.join(output, 'hub_genes.csv'))

print('Done: hubs.py')

