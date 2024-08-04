import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from scipy.stats import barnard_exact
from scipy.stats import false_discovery_control
import argparse
import os

TEST = 'fisher'

parser = argparse.ArgumentParser(description='Check all dataframes to find pathway co-occurrences')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'

t_occ_pw_files = []
n_occ_pw_files = []

print('Searching for pathways_occurrences.csv files')
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.startswith('pathways_occurrences.csv'):
            if '/tumor' in root:
                t_occ_pw_files.append(os.path.join(root, file))
            elif '/normal' in root:
                n_occ_pw_files.append(os.path.join(root, file))

def jaccard_similarity(set1, set2):
    if len(set1) == 0 or len(set2) == 0:
        return 0
    return len(set1.intersection(set2)) / len(set1.union(set2))

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

# ----- Genes overlap -----
print('Genes overlap')

def read_genesets(file):
    pw_genes = {}

    with open(file) as f:
        for line in f:
            line = line.strip().split('\t')
            # Remove empty strings and strip whitespaces
            line = [x for x in line if x]
            if line[0] not in pw_genes:
                pw_genes[line[0]] = set()
            pw_genes[line[0]].update(line[2:])

    return pw_genes

# Calculate pathways gene sets jaccard similarity

# Check if overlap file already exists
if os.path.exists('/home/lnemati/resources/reactome/reactome_overlap.csv'):
    genes_jaccard_df = pd.read_csv('/home/lnemati/resources/reactome/reactome_overlap.csv', index_col=0)
    print('Overlap file exists: /home/lnemati/resources/reactome/reactome_overlap.csv')
else:
    print('Overlap file does not exist. Creating it...')

    pw_genes = read_genesets('/home/lnemati/resources/reactome/ReactomeLowestLevelPathways.gmt')

    print('Number of pathways: ', len(pw_genes))

    print('First 3 pathways genes: ')
    for key in list(pw_genes.keys())[:3]:
        print(key)
        print(pw_genes[key])

    genes_jaccard_df = pd.DataFrame(0, index=list(pw_genes.keys()), columns=list(pw_genes.keys()))

    for i in range(len(pw_genes)):
        for j in range(i+1, len(pw_genes)):
            index1 = list(pw_genes.keys())[i]
            index2 = list(pw_genes.keys())[j]

            genes_jaccard_df.at[index1, index2] = jaccard_similarity(pw_genes[index1], pw_genes[index2])
            genes_jaccard_df.at[index2, index1] = genes_jaccard_df.at[index1, index2]
            
    # Save the overlap file
    genes_jaccard_df.to_csv('/home/lnemati/resources/reactome/reactome_overlap.csv')
    print()

all_pathways = list(genes_jaccard_df.index)

print('Number of pathways: ', len(all_pathways))
print('Number of pathways in overlap file: ', len(genes_jaccard_df))
print()

if len(all_pathways) != len(genes_jaccard_df):
    print('Error: Number of pathways in overlap file does not match the number of pathways in the current analysis')
    print('Exiting...')
    exit()

# ----- Subset -----
print('Subset')

# Subset to pairs of pathways that are linked by significant interactions
all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')

print('Reading genesets')
pw_genes = read_genesets('/home/lnemati/resources/reactome/ReactomeLowestLevelPathways.gmt')

print('Reading pathway pairs')
pairs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/reactome/reactome_connected_pathways.csv', index_col=0)

network = stack_triangle(genes_jaccard_df, 'j_genesets')
if TEST == 'fisher':
    network['log2_odds_ratio'] = 0
elif TEST == 'barnard':
    network['wald_stat'] = 0
network['pval'] = 0
network['tumor_intersection'] = 0
network['normal_intersection'] = 0
network['tumor_union'] = 0
network['normal_union'] = 0
network['keep'] = 0

all_pws = set(pairs['pathway1']).union(set(pairs['pathway2']))

# Subset to only pathways in all_pws
keep = network.index.get_level_values(0).isin(all_pws) & network.index.get_level_values(1).isin(all_pws)
network = network[keep]

# Subset to only pathways that are linked by an interaction
def are_interacting(pw1, pw2, pairs):
    result = False

    if pairs[(pairs['pathway1'] == pw1) & (pairs['pathway2'] == pw2)].shape[0] > 0:
        result = True
    if pairs[(pairs['pathway1'] == pw2) & (pairs['pathway2'] == pw1)].shape[0] > 0:
        result = True

    return result

# Only keep rows where pathway1 and pathway2 are in the same row in the pairs dataframe, order is not important
network['keep'] = network.apply(lambda x: are_interacting(x.name[0], x.name[1], pairs), axis=1)
network = network[network['keep'] == 1]

print('Number of pathways pairs connected by interactions: ', network.shape[0])

if network.shape[0] != pairs.shape[0]:
    print('Number of pathway pairs connected by interactions does not match the number of pathway pairs in the pairs dataframe')
    print('N pairs: ', pairs.shape[0]) 
    print('N network: ', network.shape[0])

# ----- Comparison -----
print('Comparison')

t_pw_modules = {}
n_pw_modules = {}

for pathway in all_pathways:
    t_pw_modules[pathway] = set()
    n_pw_modules[pathway] = set()

for path in t_occ_pw_files:
    name = path.split('/')[-2]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each pathway (idx) gets a set of all the columns where the value is > 0
    for idx in df.index:
        t_pw_modules[idx].update(df.columns[df.loc[idx] > 0])

for path in n_occ_pw_files:
    name = path.split('/')[-3]
    df = pd.read_csv(path, index_col=0)
    # add name to the columns to differentiate between different tissues
    df.columns = [name + '_' + str(x) for x in df.columns]
    # each pathway (idx) gets a set of all the columns where the value is > 0
    for idx in df.index:
        n_pw_modules[idx].update(df.columns[df.loc[idx] > 0])

print('First 3 pathways modules: ')

for key in list(t_pw_modules.keys())[:3]:
    print(key)
    print('Tumor: ', t_pw_modules[key])
    print('Normal: ', n_pw_modules[key])
    print()

min_row_col = 5

print()

# Iterate over all rows
for idx in network.index:
    pathway1 = idx[0]
    pathway2 = idx[1]
    
    network.at[idx, 'tumor_intersection'] = len(t_pw_modules[pathway1].intersection(t_pw_modules[pathway2]))
    network.at[idx, 'normal_intersection'] = len(n_pw_modules[pathway1].intersection(n_pw_modules[pathway2]))
    network.at[idx, 'tumor_union'] = len(t_pw_modules[pathway1].union(t_pw_modules[pathway2]))
    network.at[idx, 'normal_union'] = len(n_pw_modules[pathway1].union(n_pw_modules[pathway2]))
    # Fisher's exact test, make contingency table with Tumor VS Normal and Intersection VS Rest
    table = [
        [network.at[idx, 'tumor_intersection'], network.at[idx, 'tumor_union'] - network.at[idx, 'tumor_intersection']],
        [network.at[idx, 'normal_intersection'], network.at[idx, 'normal_union'] - network.at[idx, 'normal_intersection']]
    ]
    # Only keep interactions where the sum of each row and column is at least min_row_col
    if np.sum(table, axis=0)[0] < min_row_col or np.sum(table, axis=0)[1] < min_row_col:
        network.at[idx, 'keep'] = 0
    if TEST == 'fisher':
        odds_ratio, pval = fisher_exact(table)
        network.at[idx, 'log2_odds_ratio'] = np.log2(odds_ratio)
        network.at[idx, 'pval'] = pval
    elif TEST == 'barnard':
        res = barnard_exact(table)
        wald = res.statistic
        pval = res.pvalue
        network.at[idx, 'wald_stat'] = wald
        network.at[idx, 'pval'] = pval

# Convert index columns to pathway1 and pathway2
network.index = network.index.set_names(['pathway1', 'pathway2'])
network = network.reset_index()

# Create directory if it does not exist and save
if not os.path.exists(os.path.join(args.outputdir, 'pathways_network')):
    os.makedirs(os.path.join(args.outputdir, 'pathways_network'))

# Sort by statistic
if TEST == 'fisher':
    network = network.sort_values(by='log2_odds_ratio', ascending=False)
elif TEST == 'barnard':
    network = network.sort_values(by='wald_stat', ascending=False)
# Save unfiltered network
network.to_csv(os.path.join(args.outputdir, 'pathways_network', 'pathways_network_unfiltered.csv'), index=False)

# Filter
network = network[network['keep'] == 1]

# Multiple hypothesis testing correction
network = network.sort_values(by='pval')
network['pval_adj'] = false_discovery_control(network['pval'])
network = network[network['pval_adj'] < 0.05]
# Sort by statistic
if TEST == 'fisher':
    network = network.sort_values(by='log2_odds_ratio', ascending=False)
elif TEST == 'barnard':
    network = network.sort_values(by='wald_stat', ascending=False)

# Sort columns, statistic, pval_adj, pval, J_genesets, tumor_intersection, normal_intersection, tumor_union, normal_union
cols = ['pathway1', 'pathway2']
if TEST == 'fisher':
    cols += ['log2_odds_ratio']
elif TEST == 'barnard':
    cols += ['wald_stat']
cols += ['pval_adj', 'j_genesets', 'tumor_intersection', 'normal_intersection', 'tumor_union', 'normal_union']

network = network[cols]

# Save
network.to_csv(os.path.join(args.outputdir, 'pathways_network', 'pathways_network_filtered.csv'), index=False)

## ----- Tumor -----
#print('Tumor')
#
#t_pw_modules = {}
#
#for pathway in all_pathways:
#    t_pw_modules[pathway] = set()
#
#for path in t_occ_pw_files:
#    name = path.split('/')[-2]
#    df = pd.read_csv(path, index_col=0)
#    # add name to the columns to differentiate between different tissues
#    df.columns = [name + '_' + str(x) for x in df.columns]
#    # each pathway (idx) gets a set of all the columns where the value is > 0
#    for idx in df.index:
#        t_pw_modules[idx].update(df.columns[df.loc[idx] > 0])
#
#all_tumor_modules = set()
#for key in t_pw_modules.keys():
#    all_tumor_modules.update(t_pw_modules[key])
#
#print('First 3 pathways modules: ')
#for key in list(t_pw_modules.keys())[:3]:
#    print(key)
#    print(t_pw_modules[key])
#
##t_modules_jaccard_df = pd.DataFrame(0, index=all_pathways, columns=all_pathways)
#t_network = stack_triangle(genes_jaccard_df, 'j_genesets')
#t_network['j_modules'] = 0
#t_network['pathway1_modules'] = 0
#t_network['pathway2_modules'] = 0
#t_network['intersection'] = 0
#t_network['pval_adj'] = 0
#t_network['log2_odds_ratio'] = 0
#
## Iterate over all rows
#for idx in t_network.index:
#    pathway1 = idx[0]
#    pathway2 = idx[1]
#    t_network.at[idx, 'j_modules'] = jaccard_similarity(t_pw_modules[pathway1], t_pw_modules[pathway2])
#    t_network.at[idx, 'pathway1_modules'] = len(t_pw_modules[pathway1])
#    t_network.at[idx, 'pathway2_modules'] = len(t_pw_modules[pathway2])
#    t_network.at[idx, 'intersection'] = len(t_pw_modules[pathway1].intersection(t_pw_modules[pathway2]))
#
## Change index to columns pathway1 and pathway2
#t_network.index = t_network.index.set_names(['pathway1', 'pathway2'])
## Sort by adj
#t_network = t_network.sort_values(by='adj', ascending=False)
#
#print()
#
## ----- Normal -----
#print('Normal')
#
#n_pw_modules = {}
#
#for pathway in all_pathways:
#    n_pw_modules[pathway] = set()
#
#for path in n_occ_pw_files:
#    name = path.split('/')[-3]
#    df = pd.read_csv(path, index_col=0)
#    # add name to the columns to differentiate between different tissues
#    df.columns = [name + '_' + str(x) for x in df.columns]
#    # each pathway (idx) gets a set of all the columns where the value is > 0
#    for idx in df.index:
#        n_pw_modules[idx].update(df.columns[df.loc[idx] > 0])
#
#print('First 3 pathways modules: ')
#for key in list(n_pw_modules.keys())[:3]:
#    print(key)
#    print(n_pw_modules[key])
#
#n_network = stack_triangle(genes_jaccard_df, 'j_genesets')
#n_network['j_modules'] = 0
#n_network['pathway1_modules'] = 0
#n_network['pathway2_modules'] = 0
#n_network['intersection'] = 0
#n_network['adj'] = 0
#
## Iterate over all rows
#for idx in n_network.index:
#    pathway1 = idx[0]
#    pathway2 = idx[1]
#    n_network.at[idx, 'j_modules'] = jaccard_similarity(n_pw_modules[pathway1], n_pw_modules[pathway2])
#    n_network.at[idx, 'pathway1_modules'] = len(n_pw_modules[pathway1])
#    n_network.at[idx, 'pathway2_modules'] = len(n_pw_modules[pathway2])
#    n_network.at[idx, 'intersection'] = len(n_pw_modules[pathway1].intersection(n_pw_modules[pathway2]))
#    n_network.at[idx, 'adj'] = n_network.at[idx, 'j_modules'] * (1 - n_network.at[idx, 'j_genesets'])
#
## Change index to columns pathway1 and pathway2
#n_network.index = n_network.index.set_names(['pathway1', 'pathway2'])
## Sort by adj
#n_network = n_network.sort_values(by='adj', ascending=False)
#
## ----- Comparison -----
#print('Comparison')
#
## Combine the two networks
## J_genesets should be the same for both
## Add J_modules_diff column
## separate pathways1_modules and pathways2_modules into tumor and normal and add sum
## separate intersection into tumor and normal and add sum
## add adj_diff column
## make pathway1_tot and pathway2_tot columns
#
#network = t_network.copy()
#network['j_modules_tumor'] = t_network['j_modules']
#network['pathway1_modules_tumor'] = t_network['pathway1_modules']
#network['pathway2_modules_tumor'] = t_network['pathway2_modules']
#network['intersection_tumor'] = t_network['intersection']
#network['adj_tumor'] = t_network['adj']
#network.drop(columns=['j_modules', 'pathway1_modules', 'pathway2_modules', 'intersection', 'adj'], inplace=True)
#network['j_modules_normal'] = n_network['j_modules']
#network['pathway1_modules_normal'] = n_network['pathway1_modules']
#network['pathway2_modules_normal'] = n_network['pathway2_modules']
#network['intersection_normal'] = n_network['intersection']
#network['adj_normal'] = n_network['adj']
#network['j_modules_diff'] = network['j_modules_tumor'] - network['j_modules_normal']
#network['adj_diff'] = network['adj_tumor'] - network['adj_normal']
#network['intersection_tot'] = network['intersection_tumor'] + network['intersection_normal']
#network['pathway1_modules_tot'] = network['pathway1_modules_tumor'] + network['pathway1_modules_normal']
#network['pathway2_modules_tot'] = network['pathway2_modules_tumor'] + network['pathway2_modules_normal']
#
## Reorder columns
#network = network[[
#    'adj_diff',
#    'j_genesets',
#    'j_modules_diff',
#    'intersection_tot',
#    'pathway1_modules_tot',
#    'pathway2_modules_tot',
#    'j_modules_tumor',
#    'j_modules_normal',
#    'pathway1_modules_tumor',
#    'pathway1_modules_normal',
#    'pathway2_modules_tumor',
#    'pathway2_modules_normal',
#    'intersection_tumor',
#    'intersection_normal',
#    'adj_tumor',
#    'adj_normal'
#]]
#
## Sort by adj_diff
#network = network.sort_values(by='adj_diff', ascending=False)
## Save the network file
#network.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'pathways_network.csv'))
print(network)

print('Done: pathways_network.py')
