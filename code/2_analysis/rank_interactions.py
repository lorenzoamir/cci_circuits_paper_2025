import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu 
from scipy.stats import false_discovery_control
import multiprocessing
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find co-occurrences of interactors')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'

# Get number of processes from NCPUS environment variable
num_processes = int(os.getenv('NCPUS', 1))
print('Using {} processes'.format(num_processes))

t_files = []
n_files = []

all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')

print('Searching for ranked_adjacency.csv.gz files')
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.startswith('ranked_adjacency.csv'):
            if '/tumor' in root:
                t_files.append(os.path.join(root, file))
            elif '/normal' in root:
                n_files.append(os.path.join(root, file))

# ----- Comparison -----
print('Comparison')

t_values = {}
n_values = {}
interactions_genes = {}

for i, idx in enumerate(all_interactions.index):
    interaction = all_interactions.loc[idx, 'interaction']
    genes = all_interactions.loc[idx, 'interactor1':'interactor7']
    genes = list(set(genes.dropna()))
    # Skip if less than 2 genes
    if len(genes) > 1:
        interactions_genes[interaction] = genes
        t_values[interaction] = []
        n_values[interaction] = []

def triu_mean(matrix):
    return np.mean(matrix[np.triu_indices(matrix.shape[0], k=1)])

def get_value(df, genes):
    # If no genes are present, return 0
    if len(df.index.intersection(genes)) == 0:
        return 0
    # Init matrix with 0s
    matrix = pd.DataFrame(0, index=genes, columns=genes)
    # Fill the matrix with the values from the dataframe
    common_genes = list(set(genes).intersection(df.index))
    matrix.loc[common_genes, common_genes] = df.loc[common_genes, common_genes]
    return triu_mean(matrix.loc[genes, genes].values)

for path in t_files:
    df = pd.read_csv(path, index_col=0)
    for interaction in interactions_genes.keys():
        genes = interactions_genes[interaction]
        t_values[interaction].append(get_value(df, genes))

for path in n_files:
    df = pd.read_csv(path, index_col=0)
    for interaction in interactions_genes.keys():
        genes = interactions_genes[interaction]
        n_values[interaction].append(get_value(df, genes)) 

network = all_interactions.copy()
network = network[network['interaction'].isin(interactions_genes.keys())]

# Init as nan
network['U_tumor'] = np.nan
network['U_normal'] = np.nan
network['more_expressed_in'] = np.nan
network['pval'] = np.nan

print('First 3 interactions')
for key in list(t_values.keys())[:3]:
    print(key)
    print(interactions_genes[key])
    print('Tumor: ', t_values[key])
    print('Normal: ', n_values[key])
    print()

def process_chunk(chunk):
    for idx in chunk.index:
        interaction = chunk.loc[idx, 'interaction']
        U_tumor, pval = mannwhitneyu(t_values[interaction], n_values[interaction])
        U_normal = len(t_values[interaction]) * len(n_values[interaction]) - U_tumor
        chunk.at[idx, 'U_tumor'] = U_tumor
        chunk.at[idx, 'U_normal'] = U_normal
        if U_tumor > U_normal:
            chunk.at[idx, 'more_expressed_in'] = 'tumor'
        elif U_tumor < U_normal:
            chunk.at[idx, 'more_expressed_in'] = 'normal'
        chunk.at[idx, 'pval'] = pval
    return chunk

print('Splitting the dataframe into chunks')
chunks = np.array_split(network, num_processes)

print('Creating pool')
pool = multiprocessing.Pool(num_processes)

print('Processing chunks')
results = pool.map(process_chunk, chunks)

print('Closing pool')
pool.close()

print('Concatenating results')
network = pd.concat(results)
print('Network looks like this:')
print(network.head())
print()

# Create directory if it doesn't exist and save
if not os.path.exists(os.path.join(args.outputdir, 'interactions_network')):
    os.makedirs(os.path.join(args.outputdir, 'mannwhitneyu'))

# Save unfiltered network
print('Saving unfiltered network')
network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'mannwhitneyu_unfiltered.csv'), index=False)

# Filter
# Multiple hypothesis testing correction
print('Multiple hypothesis testing correction')
network = network.sort_values(by='pval')
network['pval_adj'] = false_discovery_control(network['pval'])
network = network[network['pval_adj'] < 0.05]
network = network.drop(columns=['pval'])
# Sort by statistic
network = network.sort_values(by=['more_expressed_in', 'pval_adj'], ascending=[False, True])
# Make more readable
network = network[['interaction', 'more_expressed_in', 'pval_adj']]

print(network.head())
print(network.shape)
# Save
network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'mannwhitneyu_filtered.csv'), index=False)

print('Done: rank_interactions.py')
