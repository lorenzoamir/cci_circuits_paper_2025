import pandas as pd
import numpy as np
import re
import os
import sys
import multiprocessing
import argparse
from tqdm import tqdm
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
from itertools import combinations

# Read inputfile from command line
parser = argparse.ArgumentParser(description='Find cooccurrences of mutations')

parser.add_argument('--mutationsdir', type=str, help='Input directory with mutations')
parser.add_argument('--significant', type=str, help='Input file with significant cooccurrences')

args = parser.parse_args()

# Output directory is the same as the input directory
output_dir = args.mutationsdir

# Read the input file
print('Reading input files:')
print('Significant:', args.significant)
significant = pd.read_csv(args.significant)
print('Trunk:', os.path.join(args.mutationsdir, 'Trunk', 'mutations_df.csv'))
trunk_df = os.path.join(args.mutationsdir, 'Trunk', 'mutations_df.csv')
print('Branch:', os.path.join(args.mutationsdir, 'Branch', 'mutations_df.csv'))
branch_df = os.path.join(args.mutationsdir, 'Branch', 'mutations_df.csv')
print('Private:', os.path.join(args.mutationsdir, 'Private', 'mutations_df.csv'))
private_df = os.path.join(args.mutationsdir, 'Private', 'mutations_df.csv')

# merge branch and private
trunk_df = pd.read_csv(trunk_df)
other_df = pd.concat([pd.read_csv(branch_df), pd.read_csv(private_df)])

# Remove differences between branch and private
other_df['mutation'] = 'Other'
# Strip '_Branch' and '_Private' from the patient column
other_df['module'] = other_df['module'].apply(lambda x: re.sub(r'(_Branch|_Private)', '', x))
# Strip '_Trunk' from the patient column
trunk_df['module'] = trunk_df['module'].apply(lambda x: re.sub(r'_Trunk', '', x))

all_modules = set(trunk_df['module']) | set(other_df['module'])

print('Trunk:')
print(trunk_df.head())

print('Other:')
print(other_df.head())

# Get all pairs of genes from significant cooccurrences
all_pairs = significant
all_pairs['dir_log2_odds_ratio'] = np.nan
all_pairs['dir_pval'] = np.nan

# Get number of processes from NCPUS environment variable
num_processes = int(os.getenv('NCPUS', 1))
print('Using {} processes'.format(num_processes))

def process_chunk(chunk):

    for idx, row in chunk.iterrows():
        gene1 = row['gene1']
        gene2 = row['gene2']
        
        trunk1 = set(trunk_df[trunk_df['gene'] == gene1]['module'])
        notrunk1 = all_modules - trunk1

        other2 = set(other_df[other_df['gene'] == gene2]['module'])
        noother2 = all_modules - other2

        # Make contingency table 
            
        table = np.array([
            [len(trunk1 & other2), len(trunk1 & noother2)],
            [len(notrunk1 & other2), len(notrunk1 & noother2)]
        ])
        
        # Perform Fisher's exact test
        odds_ratio, pval = fisher_exact(table, alternative='two-sided')

        chunk.at[idx, 'dir_log2_odds_ratio'] = np.log2(odds_ratio)
        chunk.at[idx, 'dir_pval'] = pval

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
all_pairs = all_pairs.sort_values(by='dir_pval')
all_pairs['dir_pval_adj'] = false_discovery_control(all_pairs['dir_pval'])

# Sort by log2 odds ratio
print('Sorting by log2 odds ratio')
all_pairs = all_pairs.sort_values('dir_log2_odds_ratio', ascending=False)

# Save unfiltered results
all_pairs.to_csv(os.path.join(output_dir, f'directionality_unfiltered.csv'), index=False)

# Filter by p-value
all_pairs = all_pairs.query('dir_pval_adj < 0.05')
all_pairs.to_csv(os.path.join(output_dir, f'directionality_filtered.csv'), index=False)

print("Done: directionality.py")
