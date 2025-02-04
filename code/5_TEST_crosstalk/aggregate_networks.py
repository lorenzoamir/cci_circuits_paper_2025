import numpy as np
import pandas as pd
import PyWGCNA
import os
import sys
import ast
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
print('Initializing network')
net = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/all_ccc_complex_pairs.csv', index_col=['complex1', 'complex2'], usecols=['complex1', 'complex2', 'all_genes'])
net['all_genes'] = net['all_genes'].apply(lambda x: ast.literal_eval(x))
net['adj'] = 0
net['same_module_frac_subtissues'] = 0
print('Initialized network:')
print(net.head())
print()

# Read list of ccc interactions and mark them in the network
ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')
# convert: string -> list -> set
ccc['all_genes'] = ccc['all_genes'].apply(lambda x: ast.literal_eval(x)).apply(set)
ccc_gene_sets = set(tuple(sorted(gene_set)) for gene_set in ccc['all_genes'])
# find the complex pairs that are actual ccc interactions
net['ccc'] = net['all_genes'].apply(lambda genes: tuple(genes) in ccc_gene_sets)
print('Number of CCC interactions: ' + str(net['ccc'].sum()))

# Read and average all networks in each subtissue directory
for subtissue in subtissue_dirs:
    print('Reading network for ' + subtissue)
    # Read subtissue network
    df = pd.read_csv(os.path.join(inputdir, subtissue, 'interactions', 'all_ccc_complex_pairs.csv'), index_col=['complex1', 'complex2'], usecols=['complex1', 'complex2', 'all_genes', 'corr'])
    # Rename corr to adj
    df = df.rename(columns={'corr': 'adj'})
    # Fill nans in adj column with 0
    df['adj'] = df['adj'].fillna(0)
    # Sum adj columns
    net['adj'] += df['adj']
    # Read wgcna file, it starts with wgcna_ and ends with .p
    wgcna = [f for f in os.listdir(os.path.join(inputdir, subtissue)) if f.startswith('wgcna_') and f.endswith('.p')][0]
    wgcna = PyWGCNA.readWGCNA(os.path.join(inputdir, subtissue, wgcna))
    # Get module memberships
    module_membership = wgcna.datExpr.var['moduleLabels']
    all_genes = set(module_membership.index)
    # If all genes are present in module_membership and are in the same module increment same_module_frac_subtissues
    for idx, net_row in net.iterrows(): 
        if all_genes.issuperset(net_row['all_genes']):
            if module_membership.loc[net_row['all_genes']].nunique() == 1:
                net.loc[idx, 'same_module_frac_subtissues'] += 1

print()

# From sum to average
net['adj'] /= len(subtissue_dirs)
net['same_module_frac_subtissues'] /= len(subtissue_dirs)

print('Aggregated network:')
print(net.head())

net.index.name = None

# Save to file
outputdir = '/home/lnemati/pathway_crosstalk/data/networks'
# Add condition to output directory
outputdir = os.path.join(outputdir, condition)
if not os.path.exists(outputdir):
    os.makedirs(outputdir)

outputfile = os.path.join(outputdir, tissue + '.csv.gz')
print('Saving to ' + outputfile)

net.to_csv(outputfile, compression='gzip')

print('Done: aggregate_networks.py')
