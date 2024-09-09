import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from itertools import combinations
import multiprocessing
import argparse
import os

# Argument parsing
parser = argparse.ArgumentParser(description='Check all dataframes to find co-occurrences of interactors')
parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')
args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

# Set directory and number of processes
parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'
num_processes = int(os.getenv('NCPUS', 1)) // 2
print(f'Using {num_processes} processes')

# Read interaction data
all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')
all_pairs = {f"{pair[0]}&{pair[1]}": [pair[0], pair[1]] for pair in combinations(all_interactions['interaction'], 2) if pair[0] != pair[1]}

# Initialize DataFrame
df = pd.DataFrame(0, index=all_pairs.keys(), columns=['N_tumor', 'N_normal'])
df['interaction1'] = df.index.map(lambda x: all_pairs[x][0])
df['interaction2'] = df.index.map(lambda x: all_pairs[x][1])

# Collect tumor and normal files
t_files, n_files = [], []
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file == 'interactions.csv':
            full_path = os.path.join(root, file)
            if '/tumor' in root:
                t_files.append(full_path)
            elif '/normal' in root:
                n_files.append(full_path)

print(f'Number of tumor files: {len(t_files)}')
print(f'Number of normal files: {len(n_files)}')

# Define co-occurrence function
def find_cooccurrences(interactions_file):
    interactions_df = pd.read_csv(interactions_file)
    interactions_df.set_index('interaction', inplace=True)
    interactions_df = interactions_df[interactions_df['same_module']]

    # Use NumPy for vectorized comparison
    cooccurrences = np.zeros(len(all_pairs))
    all_pairs_array = np.array(list(all_pairs.values()))
    valid_mask = np.isin(all_pairs_array[:, 0], interactions_df.index) & np.isin(all_pairs_array[:, 1], interactions_df.index)
    valid_pairs = all_pairs_array[valid_mask]

    if len(valid_pairs) > 0:
        module_1 = interactions_df.loc[valid_pairs[:, 0], 'module'].values
        module_2 = interactions_df.loc[valid_pairs[:, 1], 'module'].values
        matches = module_1 == module_2

        cooccurrence_indices = np.flatnonzero(valid_mask)
        cooccurrences[cooccurrence_indices] = matches.astype(int)

    return cooccurrences

# DEBUG
print('!!! DEBUG: SUBSET = 16')
SUBSET = 16
t_files = t_files[:SUBSET]
n_files = n_files[:SUBSET]

# Parallelize using imap_unordered for faster processing
with multiprocessing.Pool(num_processes) as pool:
    print('Processing tumor files...')
    t_results = list(pool.imap_unordered(find_cooccurrences, t_files))
    print('Processing normal files...')
    n_results = list(pool.imap_unordered(find_cooccurrences, n_files))

# Sum co-occurrences across all files
n_results_sum = np.sum(n_results, axis=0)
t_results_sum = np.sum(t_results, axis=0)

# Update DataFrame
df['N_tumor'] = t_results_sum
df['N_normal'] = n_results_sum

# Ensure output directory exists and save results
output_dir = os.path.join(args.outputdir, 'pairs_of_interactions')
os.makedirs(output_dir, exist_ok=True)
df.to_csv(os.path.join(output_dir, 'cooccurrences.csv'))

print('Done: pairs_of_interactions.py')

#import pandas as pd
#import numpy as np
#from scipy.stats import mannwhitneyu 
#from scipy.stats import false_discovery_control
#from itertools import combinations
#import multiprocessing
#import argparse
#import os
#
#parser = argparse.ArgumentParser(description='Check all dataframes to find co-occurrences of interactors')
#
#parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')
#
#args = parser.parse_args()
#
#if args.outputdir is None:
#    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'
#
#parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'
#
## Get number of processes from NCPUS environment variable
#num_processes = int(os.getenv('NCPUS', 1)) // 2
#print('Using {} processes'.format(num_processes))
#
#t_files = []
#n_files = []
#
#all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')
#all_pairs = {}
#
## Get all pairs
#for pair in combinations(all_interactions['interaction'], 2):
#    if pair[0] == pair[1]:
#        continue
#    pair_str = pair[0] + '&' + pair[1]
#    all_pairs[pair_str] = [pair[0], pair[1]]
#
#print('Total number of pairs: {}'.format(len(all_pairs)))
#
#print('Initializing dataframe')
#df = pd.DataFrame(0, index=all_pairs.keys(), columns=['N_tumor', 'N_normal'])
#df['interaction1'] = df.index.map(lambda x: all_pairs[x][0])
#df['interaction2'] = df.index.map(lambda x: all_pairs[x][1])
#print(df.head(3))
#print()
#
#print('Searching for interactions.csv files')
#for root, dirs, files in os.walk(parent_dir):
#    for file in files:
#        if file == 'interactions.csv':
#            if '/tumor' in root:
#                t_files.append(os.path.join(root, file))
#            elif '/normal' in root:
#                n_files.append(os.path.join(root, file))
#
#print('Number of tumor files: {}'.format(len(t_files)))
#print('Number of normal files: {}'.format(len(n_files)))
#
#def find_cooccurrences(interactions_df):
#    interactions_df = pd.read_csv(interactions_df)
#    interactions_df = interactions_df.set_index('interaction')
#
#    # Subset interactions_df to only include same module interactions
#    interactions_df = interactions_df[interactions_df['same_module']]
#
#    # Initialize the cooccurrences dictionary
#    cooccurrences = {key: 0 for key in all_pairs.keys()}
#
#    # Convert pairs to a NumPy array for efficiency
#    all_pairs_values = np.array(list(all_pairs.values()))
#
#    # Vectorized filtering: Keep only the pairs where both interactions are in the dataframe index
#    valid_mask = np.isin(all_pairs_values[:, 0], interactions_df.index) & np.isin(all_pairs_values[:, 1], interactions_df.index)
#    valid_pairs = all_pairs_values[valid_mask]
#
#    # Get the module column for both elements in the valid pairs
#    module_1 = interactions_df.loc[valid_pairs[:, 0], 'module'].values
#    module_2 = interactions_df.loc[valid_pairs[:, 1], 'module'].values
#
#    # Vectorized comparison: Set cooccurrences to 1 if the modules match, otherwise 0
#    matches = module_1 == module_2
#
#    pair_strings = [pair[0] + '&' + pair[1] for pair in valid_pairs]
#
#    # Update the cooccurrences dictionary
#    for pair_str, match in zip(pair_strings, matches):
#        cooccurrences[pair_str] = int(match)
#
#    return pd.Series(cooccurrences)
#
#print('Creating pool')
#n_pool = multiprocessing.Pool(num_processes)
#t_pool = multiprocessing.Pool(num_processes)
#
#print('Processing chunks')
#n_results = n_pool.map(find_cooccurrences, n_files)
#t_results = t_pool.map(find_cooccurrences, t_files)
#
#print('Closing pool')
#n_pool.close()
#t_pool.close()
#
#print('Results look like:')
#print(n_results[0].head(3))
#print(t_results[0].head(3))
#print()
#
#print('Checking results')
## Check that all results have the same index
#n_index = n_results[0].index
#t_index = t_results[0].index
#for i in range(1, len(n_results)):
#    assert n_index.equals(n_results[i].index)
#    assert t_index.equals(t_results[i].index)
#
#print('Concatenating results')
## Sum cooccurrences
#n_results = pd.concat(n_results, axis=1).sum(axis=1)
#t_results = pd.concat(t_results, axis=1).sum(axis=1)
#
#print('Concatenated results look like:')
#print(n_results.head(3))
#print(t_results.head(3))
#
## Save results to N_tumor and N_normal columns
#print('Saving results to dataframe')
#df['N_tumor'] = t_results
#df['N_normal'] = n_results
#print(df.head(3))
#print()
#
## Create directory if it doesn't exist and save
#os.makedirs(os.path.join(args.outputdir, 'pairs_of_interactions'), exist_ok=True)
#
## Save unfiltered network
#print('Saving') 
#df.to_csv(os.path.join(args.outputdir, 'pairs_of_interactions', 'cooccurrences.csv'))
#
#print('Done: pairs_of_interactions.py')
