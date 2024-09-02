import pandas as pd
import numpy as np
import re
import os
import multiprocessing
import argparse
from tqdm import tqdm
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
from itertools import combinations

# Read inputfile from command line
parser = argparse.ArgumentParser(description='Find cooccurrences of mutations')

parser.add_argument('--inputfile', type=str, help='Input file with mutations')

args = parser.parse_args()

# Output directory is the directory of the input file
output_dir = os.path.dirname(args.inputfile)

# Read the input file
print('Reading input file:', args.inputfile)
df = pd.read_csv(args.inputfile)
print(df.head(3))

# Iterate over gene pairs
lr_pairs = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz_all_pairs.csv', index_col=0)

# Rename columns ligand and receptor to gene1 and gene2
lr_pairs = lr_pairs.rename(columns={'ligand': 'gene1', 'receptor': 'gene2'})

# Remove rows where ligand and receptor are the same
lr_pairs = lr_pairs.query('gene1 != gene2')
print('Output directory:', output_dir)

# Subset lr_pairs to only genes that are in the df
lr_pairs = lr_pairs.query('gene1 in @df.gene and gene2 in @df.gene')

# Make a dictionary with genes as keys and modules as values
modules = df.groupby('gene')['module'].apply(set).to_dict()

print('Total number of modules:', len(set(df['module'])))
print('Number of mutated genes:', len(modules.keys()))

print('First 3 genes:')
for gene, module in list(modules.items())[:3]:
    print(f'{gene}: {module}')

all_genes = set(modules.keys())

# Subset lr_pairs to only genes that are in the modules
lr_pairs = lr_pairs.query('gene1 in @all_genes and gene2 in @all_genes')

# Generate all pairs of genes, no double counting, no self-pairs
all_pairs = list(combinations(all_genes, 2))
all_pairs = pd.DataFrame(all_pairs, columns=['gene1', 'gene2'])
all_pairs['interaction'] = 0
all_pairs['jaccard'] = 0.

# Create a MultiIndex for all_pairs
all_pairs.set_index(['gene1', 'gene2'], inplace=True)

# Create a MultiIndex for lr_pairs, including both (gene1, gene2) and (gene2, gene1)
lr_multiindex = pd.MultiIndex.from_tuples(
    list(zip(lr_pairs['gene1'], lr_pairs['gene2'])) +
    list(zip(lr_pairs['gene2'], lr_pairs['gene1']))
)

# Eeach pair is either (gene1, gene2) or (gene2, gene1), only keep the correct one
lr_multiindex = all_pairs.index[all_pairs.index.isin(lr_multiindex)]

# Update the 'interaction' column where the index matches
all_pairs.loc[lr_multiindex, 'interaction'] = 1

# Reset index if needed
all_pairs.reset_index(inplace=True)

print('all_pairs:')
print(all_pairs.head(3))

print('all_pairs shape:', all_pairs.shape)
print('Number of interacting pairs:', all_pairs['interaction'].sum())

# Get number of processes from NCPUS environment variable
num_processes = int(os.getenv('NCPUS', 1))
print('Using {} processes'.format(num_processes))

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

def process_chunk(chunk):
    # Process chunk by performing Fisher's exact test and calculating jaccard index

    for idx, row in chunk.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']
        
        # Calculate Jaccard index

        chunk.at[idx, 'jaccard'] = jaccard_similarity(set(modules[gene1]), set(modules[gene2]))

    return chunk

# Split the dataframe into chunks
print('Splitting the dataframe into chunks')
chunks = np.array_split(all_pairs, num_processes)

# Create a pool of workers
print('Creating a pool of workers')
pool = multiprocessing.Pool(num_processes)

# Compute the results
print('Computing the results')
results = pool.map(process_chunk, chunks)

# Close the pool
print('Closing the pool')
pool.close()

# Concatenate the results
print('Concatenating the results')
all_pairs = pd.concat(results)

print('Concatenated results:')
print(all_pairs)

all_pairs.to_csv(os.path.join(output_dir, 'jaccard.csv'), index=False)

print("Done: jaccard.py")
