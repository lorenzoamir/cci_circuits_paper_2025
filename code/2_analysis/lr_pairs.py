import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find co-occurrences of interactors (lr pairs)')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

# Searching for interactors_occurrences.csv files
parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'

t_occ_files = []
n_occ_files = []

all_pairs = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz_all_pairs.csv')

print('Searching for interactors_occurrences.csv files')
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.startswith('interactors_occurrences.csv'):
            if '/tumor' in root:
                t_occ_files.append(os.path.join(root, file))
            elif '/normal' in root:
                n_occ_files.append(os.path.join(root, file))

# Extend jaccard similarity to an arbitrary number of sets
def jaccard_similarity(sets):
    if len(sets) == 0:
        return 0
    if any(len(s) == 0 for s in sets):
        return 0
    intersection = set.intersection(*sets)
    union = set.union(*sets)
    return len(intersection) / len(union)

def stack_triangle(df, col):
    df = pd.DataFrame(
        df.where(
            np.tri(
                df.shape[0],
                dtype=bool,
                k=-1
            ),
            np.nan
        ).stack(dropna=True), columns=[col] 
    )
    
    return df

# ----- Comparison -----
print('Comparison')

t_modules = {}
n_modules = {}

# Init sets of modules for each gene
for i, idx in enumerate(all_pairs.index):
    genes = all_pairs.loc[idx, 'ligand':'receptor']
    genes = list(set(genes.dropna()))
    for interactor in genes:
        t_modules[interactor] = set()
        n_modules[interactor] = set()

for path in t_occ_files:
    name = path.split('/')[-2]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each gene gets a set of all the columns where the value is > 0
    for idx in df.index:
        t_modules[idx].update(df.columns[df.loc[idx] > 0])

for path in n_occ_files:
    name = path.split('/')[-2]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each pathway (idx) gets a set of all the columns where the value is > 0
    for idx in df.index:
        n_modules[idx].update(df.columns[df.loc[idx] > 0])

network = all_pairs.copy()

# Intersection is N co-occurrences
network['tumor_intersection'] = 0
network['normal_intersection'] = 0
# Union is total number of modules containing at least one gene
network['tumor_union'] = 0
network['normal_union'] = 0

network['log2_odds_ratio'] = 0
network['pval'] = 1
network['pval_adj'] = 1

network['keep'] = 1

print('First 3 interactions modules: ')
for key in list(t_modules.keys())[:3]:
    print(key)
    print('Tumor: ', t_modules[key])
    print('Normal: ', n_modules[key])
    print()


print('Calculating Jaccard similarity')
# Iterate over all rows
for idx in network.index:
    # genes are in ligand and receptor columns
    genes = network.loc[idx, 'ligand':'receptor']
    # If none of the genes are in the modules, skip
    #if not any(genes.isin(t_modules)) and not any(genes.isin(n_modules)):
    #    continue
    # Make genes into a list, make unique, remove NaNs
    genes = list(genes.dropna().unique())
    # If there is only one unique gene in the row, skip
    if len(genes) == 1:
        network.at[idx, 'keep'] = 0
        continue
    tumor_module_sets = [t_modules[gene] if gene in t_modules else set() for gene in genes]
    normal_module_sets = [n_modules[gene] if gene in n_modules else set() for gene in genes]
    network.at[idx, 'tumor_intersection'] = len(set.intersection(*tumor_module_sets))
    network.at[idx, 'normal_intersection'] = len(set.intersection(*normal_module_sets))
    network.at[idx, 'tumor_union'] = len(set.union(*tumor_module_sets))
    network.at[idx, 'normal_union'] = len(set.union(*normal_module_sets))
    # Make contingency table with Tumor VS Normal and Intersection VS Rest
    table = [
        [network.at[idx, 'tumor_intersection'], network.at[idx, 'tumor_union'] - network.at[idx, 'tumor_intersection']],
        [network.at[idx, 'normal_intersection'], network.at[idx, 'normal_union'] - network.at[idx, 'normal_intersection']]
    ]
    odds_ratio, pval = fisher_exact(table)
    network.at[idx, 'log2_odds_ratio'] = np.log2(odds_ratio)
    network.at[idx, 'pval'] = pval

# Create directory if it doesn't exist and save
if not os.path.exists(os.path.join(args.outputdir, 'lr_pairs')):
    os.makedirs(os.path.join(args.outputdir, 'lr_pairs'))

# Sort by statistic
network = network.sort_values(by='log2_odds_ratio', ascending=False)
# Save unfiltered network
network.to_csv(os.path.join(args.outputdir, 'lr_pairs', 'lr_pairs_unfiltered.csv'), index=False)

# Filter
network = network[network['keep'] == 1]
# Multiple hypothesis testing correction
network = network.sort_values(by='pval')
network['pval_adj'] = false_discovery_control(network['pval'])
network = network[network['pval_adj'] < 0.05]
# Sort by statistic
network = network.sort_values(by='log2_odds_ratio', ascending=False)

# Sort columns, statistic, pval_adj, pval, J_genesets, tumor_intersection, normal_intersection, tumor_union, normal_union
cols = ['interaction']
cols += ['log2_odds_ratio']
cols += ['pval_adj', 'tumor_intersection', 'normal_intersection', 'tumor_union', 'normal_union']

print(network.head())
print(network.shape)
# Save
network.to_csv(os.path.join(args.outputdir, 'lr_pairs', 'lr_pairs_filtered.csv'), index=False)

print('Done: lr_pairs.py')
