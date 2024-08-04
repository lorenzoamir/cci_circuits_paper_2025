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
# Subset df to only genes that are in the lr_pairs
lr_genes = set(lr_pairs['gene1']).union(set(lr_pairs['gene2']))

# Make a dictionary with genes as keys and modules as values
modules = df.groupby('gene')['module'].apply(set).to_dict()
all_modules = set(df['module'])

print('Total number of modules:', len(all_modules))
print('Number of mutated genes:', len(modules.keys()))

# Subset to only lr_genes
modules = {gene: modules[gene] for gene in lr_genes}
print('Number of mutated genes in lr_pairs:', len(modules.keys()))

print('First 3 genes:')
for gene, module in list(modules.items())[:3]:
    print(f'{gene}: {module}')

# Subset to only genes with at least min_modules modules
min_modules = 0
print(f'Subsetting to genes mutated in at least {min_modules} patients')
all_genes = set([gene for gene in modules.keys() if len(modules[gene]) >= min_modules])
print(f'Number of genes: {len(all_genes)}')
modules = {gene: modules[gene] for gene in all_genes}
# Subset to only all_genes in lr_pairs
lr_pairs = lr_pairs.query('gene1 in @all_genes and gene2 in @all_genes')

# Get all pairs from lr_pairs
all_pairs = list(zip(lr_pairs['gene1'], lr_pairs['gene2']))

# Make all_pairs a DataFrame
all_pairs = pd.DataFrame(all_pairs, columns=['gene1', 'gene2'])

print(f'Total number of pairs: {len(all_pairs)}')

# Add jaccard, log2_odds_ratio and p_value columns
all_pairs['jaccard'] = np.nan
all_pairs['log2_odds_ratio'] = np.nan
all_pairs['pval'] = np.nan

# Add empty sets for genes that are not in the modules
#for gene in lr_genes.difference(modules.keys()):
#    modules[gene] = set()

# Get number of processes from NCPUS environment variable
num_processes = int(os.getenv('NCPUS', 1))
print('Using {} processes'.format(num_processes))

def process_chunk(chunk):
    # Process chunk by performing Fisher's exact test and calculating jaccard index

    for idx, row in chunk.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']

        # Make contingency table, one axis is module contains gene1 or not,
        # the other is module contains gene2 or not

        # Contingency table
        gene1_with = set(modules[gene1])
        gene1_without = all_modules - gene1_with
        gene2_with = set(modules[gene2])
        gene2_without = all_modules - gene2_with

        # Make contingency table
        table = np.array([
            [len(gene1_with & gene2_with), len(gene1_with & gene2_without)],
            [len(gene1_without & gene2_with), len(gene1_without & gene2_without)]
        ])
        
        # Calculate Jaccard index
        jaccard = len(gene1_with & gene2_with) / len(gene1_with | gene2_with)
        # Perform Fisher's exact test
        odds_ratio, pval = fisher_exact(table, alternative='two-sided')

        chunk.at[idx, 'jaccard'] = jaccard
        chunk.at[idx, 'log2_odds_ratio'] = np.log2(odds_ratio)
        chunk.at[idx, 'pval'] = pval

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

# Adjust p-values
print('Adjusting p-values')
all_pairs = all_pairs.sort_values(by='pval')
all_pairs['pval_adj'] = false_discovery_control(all_pairs['pval'])

# Sort by log2 odds ratio
print('Sorting by log2 odds ratio')
all_pairs = all_pairs.sort_values('log2_odds_ratio', ascending=False)

# Save unfiltered results
all_pairs.to_csv(os.path.join(output_dir, f'cooccurrences_unfiltered.csv'), index=False)

# Filter by p-value
all_pairs = all_pairs.query('pval_adj < 0.05')
all_pairs.to_csv(os.path.join(output_dir, f'cooccurrences_filtered.csv'), index=False)

print("Done: cooccurrence.py")
