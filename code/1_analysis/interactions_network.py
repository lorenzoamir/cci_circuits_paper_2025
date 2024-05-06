import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find co-occurrences of interactors')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

# Searching for interactors_occurrences.csv files
parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'

t_occ_files = []
n_occ_files = []

all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')

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

for i, idx in enumerate(all_interactions.index):
    interaction = all_interactions.loc[idx, 'interaction']
    genes = all_interactions.loc[idx, 'interactor1':'interactor7']
    genes = list(set(genes.dropna()))
    # DEBUG:
    if i < 5:
        print(interaction)
        print(genes)
        print()
    for interactor in genes:
        t_modules[interactor] = set()
        n_modules[interactor] = set()

for path in t_occ_files:
    name = path.split('/')[-2]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each pathway (idx) gets a set of all the columns where the value is > 0
    for idx in df.index:
        t_modules[idx].update(df.columns[df.loc[idx] > 0])

for path in n_occ_files:
    name = path.split('/')[-3]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each pathway (idx) gets a set of all the columns where the value is > 0
    for idx in df.index:
        n_modules[idx].update(df.columns[df.loc[idx] > 0])

network = all_interactions.copy()
network['log2_odds_ratio'] = 0
network['pval'] = 0
network['tumor_intersection'] = 0
network['normal_intersection'] = 0
network['tumor_union'] = 0
network['normal_union'] = 0
network['keep'] = 1 

print('First 3 interactions modules: ')
for key in list(t_modules.keys())[:3]:
    print(key)
    print('Tumor: ', t_modules[key])
    print('Normal: ', n_modules[key])
    print()

min_intersection = 15
min_union = 15

print()

print('Iterating over all interactions')
# Iterate over all rows
for idx in network.index:
    genes = network.loc[idx, 'interactor1':'interactor7']
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
    # Only keep interactions with at least min_intersection and min_union
    if network.at[idx, 'tumor_intersection'] + network.at[idx, 'normal_intersection'] < min_intersection:
        network.at[idx, 'keep'] = 0
        continue
    if network.at[idx, 'tumor_union'] < min_union or network.at[idx, 'normal_union'] < min_union:
        network.at[idx, 'keep'] = 0
        continue
    # Fisher's exact test, make contingency table with Tumor VS Normal and Intersection VS Rest
    table = [
        [network.at[idx, 'tumor_intersection'], network.at[idx, 'tumor_union'] - network.at[idx, 'tumor_intersection']],
        [network.at[idx, 'normal_intersection'], network.at[idx, 'normal_union'] - network.at[idx, 'normal_intersection']]
    ]
    odds_ratio, pval = fisher_exact(table)
    network.at[idx, 'log2_odds_ratio'] = np.log2(odds_ratio)
    network.at[idx, 'pval'] = pval

# Create directory if it doesn't exist and save
if not os.path.exists(os.path.join(args.outputdir, 'interactions_network')):
    os.makedirs(os.path.join(args.outputdir, 'interactions_network'))

# Sort by log2_odds_ratio
network = network.sort_values(by='log2_odds_ratio', ascending=False)
# Save unfiltered network
network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'interactions_network_unfiltered.csv'), index=False)

# Filter
network = network[network['keep'] == 1]
# Multiple hypothesis testing correction
network = network.sort_values(by='pval')
network['pval_adj'] = false_discovery_control(network['pval'])
# Sort 
network = network.sort_values(by='log2_odds_ratio', ascending=False)
# Sort columns, log2_odds_ratio, pval_adj, pval, J_genesets, tumor_intersection, normal_intersection, tumor_union, normal_union
network = network[[
    'interaction',
    'log2_odds_ratio', 'pval_adj',
    'tumor_intersection', 'normal_intersection',
    'tumor_union', 'normal_union'
]]

# Filter and save
network = network[network['pval_adj'] < 0.05]
network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'interactions_network_filtered.csv'), index=False)

## ----- Tumor -----
#print('Tumor')
#
#t_modules = {}
#
#for path in t_occ_files:
#    name = path.split('/')[-2]
#    df = pd.read_csv(path, index_col=0)
#    # add name to the columns to differentiate between different tissues
#    df.columns = [name + '_' + str(x) for x in df.columns]
#    # each pathway (idx) gets a set of all the columns where the value is > 0
#    for idx in df.index:
#        if idx not in t_modules:
#            t_modules[idx] = set()
#        t_modules[idx].update(df.columns[df.loc[idx] > 0])
#
#print('First 3 genes modules: ')
#for key in list(t_modules.keys())[:3]:
#    print(key)
#    print(t_modules[key])
#
##t_modules_jaccard_df = pd.DataFrame(0, index=all_pathways, columns=all_pathways)
#t_network = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')
#t_network['j_modules'] = 0
#t_network['intersection'] = 0
#t_network['union'] = 0
#for i in range(1, 8):
#    t_network['interactor' + str(i) + '_modules'] = 0
#t_network['keep'] = 1
#
## Iterate over all rows
#for idx in t_network.index:
#    # Genes are the genes in cols interactor1, interactor2, ..., interactor7
#    genes = t_network.loc[idx, 'interactor1':'interactor7']
#    # If none of the genes are in the modules, skip
#    if not any(genes.isin(t_modules)):
#        continue
#    # Make genes into a list, remove NaNs
#    genes = list(set(genes.dropna()))
#    # If there is only one unique gene in the row, mark row to be removed and skip
#    if len(genes) == 1:
#        t_network.at[idx, 'keep'] = 0
#        continue
#    module_sets = [t_modules[gene] if gene in t_modules else set() for gene in genes]
#    t_network.at[idx, 'j_modules'] = jaccard_similarity(module_sets)
#    t_network.at[idx, 'intersection'] = len(set.intersection(*module_sets))
#    t_network.at[idx, 'union'] = len(set.union(*module_sets))
#    for i, gene in enumerate(genes):
#        if gene in t_modules:
#            t_network.at[idx, 'interactor' + str(i+1) + '_modules'] = len(t_modules[gene])
#
## Remove rows with only one unique gene
#t_network = t_network[t_network['keep'] == 1]
#t_network.drop(columns=['keep'], inplace=True)
## Sort by j_modules and intersection
#t_network = t_network.sort_values(by=['j_modules', 'intersection'], ascending=False)
## Save the network file
#os.makedirs(os.path.join(args.outputdir, 'interactions_network'), exist_ok=True)
#t_network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'tumor_interactions_network.csv'), index=False)
#
#print()
#
## ----- Normal -----
#print('Normal')
#
#n_modules = {}
#
#for path in n_occ_files:
#    name = path.split('/')[-3]
#    df = pd.read_csv(path, index_col=0)
#    # add name to the columns to differentiate between different tissues
#    df.columns = [name + '_' + str(x) for x in df.columns]
#    # each pathway (idx) gets a set of all the columns where the value is > 0
#    for idx in df.index:
#        if idx not in n_modules:
#            n_modules[idx] = set()
#        n_modules[idx].update(df.columns[df.loc[idx] > 0])
#
#print('First 3 genes modules: ')
#for key in list(n_modules.keys())[:3]:
#    print(key)
#    print(n_modules[key])
#
#n_network = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')
#n_network['j_modules'] = 0
#n_network['intersection'] = 0
#n_network['union'] = 0
#for i in range(1, 8):
#    n_network['interactor' + str(i) + '_modules'] = 0
#n_network['keep'] = 1
#
## Iterate over all rows
#for idx in n_network.index:
#    # Genes are the genes in cols interactor1, interactor2, ..., interactor7
#    genes = n_network.loc[idx, 'interactor1':'interactor7']
#    # If none of the genes are in the modules, skip
#    if not any(genes.isin(n_modules)):
#        continue
#    # Make genes into a list, remove NaNs
#    genes = list(set(genes.dropna()))
#    # If there is only one unique gene in the row, mark row to be removed and skip
#    if len(genes) == 1:
#        n_network.at[idx, 'keep'] = 0
#        continue
#    module_sets = [n_modules[gene] if gene in n_modules else set() for gene in genes]
#    n_network.at[idx, 'j_modules'] = jaccard_similarity(module_sets)
#    n_network.at[idx, 'intersection'] = len(set.intersection(*module_sets))
#    n_network.at[idx, 'union'] = len(set.union(*module_sets))
#    for i, gene in enumerate(genes):
#        if gene in n_modules:
#            n_network.at[idx, 'interactor' + str(i+1) + '_modules'] = len(n_modules[gene])
#
## Remove rows with only one unique gene
#n_network = n_network[n_network['keep'] == 1]
#n_network.drop(columns=['keep'], inplace=True)
## Sort by j_modules and intersection
#n_network = n_network.sort_values(by=['j_modules', 'intersection'], ascending=False)
## Save the network file
#n_network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'normal_interactions_network.csv'), index=False)
#
## ----- Comparison -----
#print('Comparison')
#
## Combine the two networks
## Add J_modules_diff column
## separate interactors_modules into tumor and normal sum
## separate intersection into tumor and normal and add sum
## make interactorN_tot columns
#t_network = t_network.set_index('interaction')
#n_network = n_network.set_index('interaction')
#
#network = t_network.copy()
#network['j_modules_tumor'] = t_network['j_modules']
#network['j_modules_normal'] = n_network.loc[t_network.index, 'j_modules']
#network['intersection_tumor'] = t_network['intersection']
#network['intersection_normal'] = n_network.loc[t_network.index, 'intersection']
#network['union_tumor'] = t_network['union']
#network['union_normal'] = n_network.loc[t_network.index, 'union']
#network['j_modules_diff'] = network['j_modules_tumor'] - network['j_modules_normal']
#network['intersection_tot'] = network['intersection_tumor'] + network['intersection_normal']
#network['union_tot'] = network['union_tumor'] + network['union_normal']
#
#
## Reorder columns
#network = network[[
#    'j_modules_diff', 'j_modules_tumor', 'j_modules_normal',
#    'intersection_tot', 'intersection_tumor', 'intersection_normal',
#    'union_tot', 'union_tumor', 'union_normal',
#]]
#
## Sort by j_modules_diff
#network = network.sort_values(by='j_modules_diff', ascending=False)
## Save the network file
#network.to_csv(os.path.join(args.outputdir, 'interactions_network', 'interactions_network.csv'), index=True)

print('Done: pathways_cooccurrences.py')
