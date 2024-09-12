import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from itertools import combinations
import multiprocessing as mp
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

# ----- Process files -----

# Main function for multiprocessing
def parallel_processing(all_paths, num_processes=1):
    
    all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv', index_col='interaction')
    cooc_result = pd.DataFrame(0, index=all_interactions.index, columns=all_interactions.index)
    occ_result = pd.DataFrame(0, index=all_interactions.index, columns=all_interactions.index)

    all_interactions['module'] = np.nan
    all_interactions = all_interactions['module']

    # Create a pool of workers
    with mp.Pool(num_processes) as pool:
        # Map the process_chunk function to each chunk
        cooc_occ_tuples = pool.map(coocs_occs, [(path, all_interactions) for path in all_paths])

    # Combine results by summing all coocs matrices to results
    for cooc, occ in cooc_occ_tuples:
        cooc_result = cooc_result.add(cooc.astype(int), fill_value=np.nan)
        occ_result = occ_result.add(occ.astype(int), fill_value=np.nan)

    # Calculate Jaccard index, its just the number of co-occurrences divided by the number of total occurrences
    jacc_result = cooc_result.div(occ_result, fill_value=0)
    jacc_result = jacc_result.fillna(0) # If no occurrences, the Jaccard index is 0
    
    # Fill diagonal with 0
    np.fill_diagonal(cooc_result.values, 0)
    np.fill_diagonal(jacc_result.values, 0)

    cooc_result.index.name = None
    occ_result.index.name = None
    jacc_result.index.name = None

    return cooc_result, occ_result, jacc_result

def coocs_occs(args):
    # Calculate the co-occurrences matrix, the total number of occurrences, and the Jaccard index

    path, all_interactions = args
    # Read the file
    modules = pd.read_csv(path, index_col='interaction')['module']
    # Copy all_interactions to avoid modifying the original
    all_interactions = all_interactions.copy()
    all_interactions.loc[modules.index.intersection(all_interactions.index)] = modules

    # Co-occurrences is the number of times the interactions appear together in the same module
    coocs = pd.DataFrame(np.equal.outer(all_interactions.values, all_interactions.values), index=all_interactions.index, columns=all_interactions.index)
    
    # Total occs is the number of times the interactions appear in a module, so it can be 0 (None), 1 (one of them), or 2 (both)
    valid = (~all_interactions.isna()).astype(int).values # Get interactions that actually have a module
    total_occs = pd.DataFrame(np.add.outer(valid, valid), index=all_interactions.index, columns=all_interactions.index)
    # We are just counting the number of times each interaction appears in a module, but if they are in the same module, we have to subtract 1
    total_occs = total_occs - coocs

    # Return matrices
    return (coocs, total_occs)

# Process tumor and normal files separately
t_coocs, t_occs, t_jacc = parallel_processing(t_files, num_processes)
n_coocs, n_occs, n_jacc = parallel_processing(n_files, num_processes)

# Ensure output directory exists and save results
output_dir = os.path.join(args.outputdir, 'pairs_of_interactions')
os.makedirs(output_dir, exist_ok=True)

t_coocs.to_csv(os.path.join(output_dir, 'tumor_coocs.csv'))
n_coocs.to_csv(os.path.join(output_dir, 'normal_coocs.csv'))

t_occs.to_csv(os.path.join(output_dir, 'tumor_occs.csv'))
n_occs.to_csv(os.path.join(output_dir, 'normal_occs.csv'))

t_jacc.to_csv(os.path.join(output_dir, 'tumor_jaccard.csv'))
n_jacc.to_csv(os.path.join(output_dir, 'normal_jaccard.csv'))

print('Done: pairs_of_interactions.py')
