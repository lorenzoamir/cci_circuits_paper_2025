import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
from scipy.stats import mannwhitneyu, false_discovery_control
from multiprocessing import Pool
import ast
import sys

tcolor     = '#ab3502'
ncolor     = '#00728e'
graycolor  = '#4D4E4F'
graycolor2 = '#C8CAD4'

data_dir='/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/'

# Find all subdirectories of datadir, those are the tissues
tissue_dirs = [os.path.join(data_dir, x) for x in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, x))]
print('Tissues:', len(tissue_dirs))

# Tissue directories have tumor and normal subdirectories, each has subdirectories with subtissues
tumor_dirs = {}
normal_dirs = {}

for tissue_dir in tissue_dirs:
    # Get tissue name from path
    tissue_name = tissue_dir.split("/")[-1]
    
    # Find normal and tumor subdirectories
    tumor_dirs[tissue_name] = [os.path.join(tissue_dir, 'tumor', x) for x in os.listdir(os.path.join(tissue_dir, 'tumor'))]
    normal_dirs[tissue_name] = [os.path.join(tissue_dir, 'normal', x) for x in os.listdir(os.path.join(tissue_dir, 'normal'))]
    
    # Exclude empty subtissues dirs that don't have any files inside
    tumor_dirs[tissue_name] = [x for x in tumor_dirs[tissue_name] if os.listdir(x)]
    normal_dirs[tissue_name] = [x for x in normal_dirs[tissue_name] if os.listdir(x)]

# Make a dataframe with tissue, subtissue, condition
tissues_df = pd.DataFrame(columns=['tissue', 'subtissue', 'condition'])

# Create dataframes for tumor and normal conditions
tumor_data = [{'tissue': tissue_name, 'subtissue': subtissue, 'condition': 'tumor'}
              for tissue_name in tumor_dirs for subtissue in tumor_dirs[tissue_name]]
normal_data = [{'tissue': tissue_name, 'subtissue': subtissue, 'condition': 'normal'}
               for tissue_name in normal_dirs for subtissue in normal_dirs[tissue_name]]

# Concatenate the dataframes
tissues_df = pd.concat([tissues_df, pd.DataFrame(tumor_data), pd.DataFrame(normal_data)], ignore_index=True)


# Store directories in all_tissues
tissues_df['path'] = tissues_df.subtissue
tissues_df['subtissue'] = tissues_df.subtissue.apply(lambda x: x.split('/')[-1])

# Print each path and whether it exists
for path in tissues_df['path']:
    print(path, os.path.exists(path))

print('Tissues dataframe:', tissues_df.shape)
print(tissues_df.head())

# Save tissues_df to a csv file
tissues_df.to_csv('/home/lnemati/pathway_crosstalk/data/tissues.csv', index=False)

# Make new normal and tumor dirs with new tissue names as keys and paths as values
tumor_dirs = {}
normal_dirs = {}

for i, row in tissues_df.iterrows():
    if row.condition == 'tumor':
        if row.tissue not in tumor_dirs:
            tumor_dirs[row.tissue] = []
        tumor_dirs[row.tissue].append(row.path)
    else:
        if row.tissue not in normal_dirs:
            normal_dirs[row.tissue] = []
        normal_dirs[row.tissue].append(row.path)

print('Tumor dirs:', tumor_dirs)
print('Normal dirs:', normal_dirs)

#filename = 'all_ccc_gene_pairs.csv'
filename = 'all_ccc_complex_pairs.csv'
#metrics = ['adj', 'kme_corr', 'corr']
metrics = ['adj']

print('filename:', filename)

all_interactions = pd.read_csv(os.path.join('/home/lnemati/pathway_crosstalk/data/interactions', filename), index_col='interaction')

# Init columns for each metric
for metric in metrics:
    all_interactions[f'{metric}_normal'] = 0
    all_interactions[f'{metric}_tumor'] = 0

normal_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}
tumor_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}
missing_genes = pd.DataFrame(index=all_interactions.index)

def process_tissue(row, metric, filename, all_interactions):
    tissue = row.tissue
    path = row.path
    subtissue = row.subtissue
    condition = row.condition

    # Make a copy of all_interactions for this tissue
    tissue_interactions = all_interactions.copy()
    tissue_interactions[metric] = 0

    # Read the interactions table for this subtissue
    df = pd.read_csv(os.path.join(path, 'interactions', filename), index_col='interaction')
    df = df.fillna(0)

    # Add the metric values to the tissue_interactions dataframe
    tissue_interactions[metric] = df[metric]

    # Prepare result dictionary to update main data structures in parent process
    result = {
        'tissue': tissue,
        'subtissue': subtissue,
        'condition': condition,
        'metric_values': tissue_interactions[metric],
        'missing_genes': df['missing_genes']
    }
    return result

def parallel_process(metric, tissues_df, filename, all_interactions):
    ncpus = os.environ.get('NCPUS')
    with Pool(int(ncpus)) as pool:
        results = pool.starmap(
            process_tissue,
            [(row, metric, filename, all_interactions) for _, row in tissues_df.iterrows()]
        )

    # Aggregate results into normal_vals, tumor_vals, and missing_genes
    for result in results:
        subtissue = result['subtissue']
        condition = result['condition']
        metric_values = result['metric_values']
        missing_genes[subtissue] = result['missing_genes']

        if condition == 'normal':
            normal_vals[metric][subtissue] = np.nan
            normal_vals[metric].loc[metric_values.index, subtissue] = metric_values
        elif condition == 'tumor':
            tumor_vals[metric][subtissue] = np.nan
            tumor_vals[metric].loc[metric_values.index, subtissue] = metric_values

    # Save the values per tissue to a csv file
    outdir = os.path.join('/home/lnemati/pathway_crosstalk/results/crosstalk', filename.replace('.csv', ''), metric)
    os.makedirs(outdir, exist_ok=True)
    print('Saving: ', f'{outdir}/normal_values.csv')
    normal_vals[metric].to_csv(f'{outdir}/normal_values.csv')
    print('Saving: ', f'{outdir}/tumor_values.csv')
    tumor_vals[metric].to_csv(f'{outdir}/tumor_values.csv')
    print('Saving: ', f'{outdir}/missing_genes.csv')
    missing_genes.to_csv(f'{outdir}/missing_genes.csv')

# Main loop to parallelize over metrics
for metric in metrics:
    parallel_process(metric, tissues_df, filename, all_interactions)

# Merge results together
tissues_df = pd.read_csv('/home/lnemati/pathway_crosstalk/data/tissues.csv')
tissues_df.loc['thyroid_carcinoma', 'tissue'] = 'thyroid'
tissues_df.index = tissues_df.subtissue

paths2sub = dict(tissues_df.path)
paths2sub = {v: k for k, v in paths2sub.items()}
sub2tissue = dict(tissues_df.tissue)
sub2condition = dict(tissues_df.condition)

resource = 'all_ccc_complex_pairs'
metric = 'adj'
alternative='two-sided'

outdir = f'/home/lnemati/pathway_crosstalk/results/crosstalk/{resource}/{metric}'

# Read normal and tumor values
normal_vals = pd.read_csv(os.path.join(outdir, 'normal_values.csv'), index_col=0)
tumor_vals = pd.read_csv(os.path.join(outdir, 'tumor_values.csv'), index_col=0)

# Group subtissues by major tissue
normal_vals = normal_vals.T.groupby(lambda x: sub2tissue[x]).mean().T # group normal subtissues by major tissue
tumor_vals = tumor_vals.T.groupby(lambda x: sub2tissue[x]).mean().T # group tumor subtissues by major tissue

# Get pairs and find those sharing genes
all_pairs = normal_vals.index

# get series for complex1 and complex2
print('Finding pairs with shared genes')
c1 = pd.Series(all_pairs.str.split('+').str[0], index=all_pairs) # Get complexes from pairs
c1 = c1.apply(lambda x: x.split('_')) # Split genes in each complex
c2 = pd.Series(all_pairs.str.split('+').str[1], index=all_pairs)
c2 = c2.apply(lambda x: x.split('_'))

common = c1.combine(c2, lambda x, y: set(x) & set(y)) # Find common genes between complexes
common = common.apply(lambda x: len(x) > 0) # Check if there are common genes
print('Number of pairs with shared genes:', common.sum())

# Get standard deviation (used for filtering)
std = normal_vals.merge(tumor_vals, left_index=True, right_index=True).std(axis=1)

# Don't test pairs with all 0 values or pairs that share genes
min_std = 0

print('Number of pairs with low variance:', (std <= min_std).sum())

# Get test results
results = mannwhitneyu(tumor_vals.T, normal_vals.T, alternative=alternative)
pvals = results.pvalue
auroc = results.statistic / (tumor_vals.shape[1] * normal_vals.shape[1])

pvals = pd.Series(pvals, all_pairs)
auroc = pd.Series(auroc, all_pairs)

# Make network
network = pd.DataFrame(auroc, index=all_pairs, columns=['auroc'])
network['pval'] = pvals

# Add standard deviation
network['std'] = normal_vals.merge(tumor_vals, left_index=True, right_index=True).std(axis=1)
# Add normal and tumor median
network['n_median'] = normal_vals.median(axis=1)
network['t_median'] = tumor_vals.median(axis=1)

# Correct pvalues, don't consider those that have genes in common
common_interactions = common[common].index
no_std = network[std <= min_std].index

no_test = network.index.intersection(common_interactions.union(no_std))
yes_test = network.index.difference(no_test)

# Set pval of common and no_std interactions to None, correct the rest
network.loc[no_test, 'pval'] = None
network.loc[no_test, 'pval_adj'] = None
network.loc[yes_test, 'pval_adj'] = false_discovery_control(network.loc[yes_test, 'pval'])
#network['pval_adj'] = false_discovery_control(pvals)
network.head()

network['all_genes'] = network.index.str.replace('+', '_').str.split('_')
network['all_genes'] = network['all_genes'].apply(lambda x: tuple(sorted(x)))

network.insert(0, 'complex1', network.index.str.split('+', expand=True).get_level_values(0))
network.insert(1, 'complex2', network.index.str.split('+', expand=True).get_level_values(1))

# Significant pairs
significant_pairs = network.query('pval_adj < 0.05').index
print('Number of significant pairs:', len(significant_pairs))

# Add ccc information
ccc = pd.read_csv('/home/lnemati/pathway_crosstalk/data/interactions/ccc.csv')

# convert: string -> list -> set
ccc['all_genes'] = ccc['all_genes'].apply(lambda x: ast.literal_eval(x)).apply(set)
ccc_gene_sets = set(tuple(sorted(gene_set)) for gene_set in ccc['all_genes'])

# find the complex pairs that are actual ccc interactions
network['ccc'] = network['all_genes'].apply(lambda genes: tuple(genes) in ccc_gene_sets)

print('Number of CCC pairs (significant and not):', network.query('ccc').shape[0])

print('probability of ccc interactions being significant')
print(network.query('ccc and pval_adj < 0.05').shape[0] / network.query('ccc').shape[0])

print('probability of other complex pairs being significant')
print(network.query('(not ccc) and pval_adj < 0.05').shape[0] / network.query('not ccc').shape[0])

# Save network

network = network.sort_values(by='pval')
network.to_csv(os.path.join(outdir, 'network_unfiltered.csv.gz'), index=False)
network.query('pval_adj < 0.05').to_csv(os.path.join(outdir, 'network_filtered.csv'), index=False)
