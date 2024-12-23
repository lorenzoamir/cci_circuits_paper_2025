import numpy as np
import pandas as pd
import os
import sys
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Aggregate subtissue networks into tissue networks')
parser.add_argument('--inputdir', type=str, help='Input directory containing subtissue networks')

args = parser.parse_args()

inputdir = args.inputdir
print('Input directory: ' + inputdir)

# Read tissue and condition from input directory
tissue = os.path.basename(os.path.dirname(inputdir))
condition = os.path.basename(inputdir)
print('Tissue: ' + tissue)
print('Condition: ' + condition)

# Get all directories in input directory
subtissue_dirs = [f for f in os.listdir(inputdir) if os.path.isdir(os.path.join(inputdir, f))]
print('Subtissues: ' + ', '.join(subtissue_dirs))

# Initialize empty adjacency matrix
print('Initializing adjacency matrix')
adj = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/all_ccc_complex_pairs.csv', index_col=['complex1', 'complex2'], usecols=['complex1', 'complex2'])
# Make symmetric by concatenating with a df with complex1 and complex2 (multi index) swapped
adj = pd.concat([adj, adj.reset_index().rename(columns={'complex1': 'complex2', 'complex2': 'complex1'}).set_index(['complex1', 'complex2'])])
adj['adj'] = 0
adj = adj['adj'].unstack(fill_value=0).fillna(0)
print('Initialized adjacency matrix:')
print(adj.head())
print()

# Read and average all networks in each subtissue directory
for subtissue in subtissue_dirs:
    print('Reading network for ' + subtissue)
    # Read network and convert to adjacency matrix
    df = pd.read_csv(os.path.join(inputdir, subtissue, 'interactions', 'all_ccc_complex_pairs.csv'), index_col=['complex1', 'complex2'], usecols=['complex1', 'complex2', 'adj'])
    # Make symmetric
    df_swp = df.reset_index().rename(columns={'complex1': 'complex2', 'complex2': 'complex1'}).set_index(['complex1', 'complex2'])
    print(df.head())
    print(df_swp.head())
    df = pd.concat([df, df_swp])
    df = df['adj'].unstack(fill_value=0).fillna(0)
    print(df.head())
    # Add to aggregate adjacency matrix
    adj = adj.add(df, fill_value=0)

print()

# From sum to average
adj = adj / len(subtissue_dirs)

print('Aggregated adjacency matrix:')
print(adj.head())

# Save to file
outputdir = '/home/lnemati/pathway_crosstalk/data/networks'
# Add condition to output directory
outputdir = os.path.join(outputdir, condition)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

outputfile = os.path.join(outputdir, tissue + '.csv.gz')
print('Saving to ' + outputfile)

adj.to_csv(outputfile, compression='gzip')

print('Done: aggregate_networks.py')
