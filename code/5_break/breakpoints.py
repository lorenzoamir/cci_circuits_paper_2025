import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
from multiprocessing import Pool
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
tissues_df = tissues_df.append([{'tissue': tissue_name, 'subtissue': subtissue, 'condition': 'tumor'} for tissue_name in tumor_dirs for subtissue in tumor_dirs[tissue_name]])
tissues_df = tissues_df.append([{'tissue': tissue_name, 'subtissue': subtissue, 'condition': 'normal'} for tissue_name in normal_dirs for subtissue in normal_dirs[tissue_name]])

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

# Find tissues that have both tumor and normal conditions
#tissues = tissues_df.groupby('tissue').condition.nunique()
#tissues = tissues[tissues == 2].index
#print('Tissues with both tumor and normal:', len(tissues))

# Remove normal testis subtissues
# tissues_df = tissues_df[~((tissues_df.tissue == 'testis') & (tissues_df.condition == 'normal'))]
# tissues = [x for x in tissues if x != 'testis']

# Join thyroid_gland and thyroid subtissues into the same tissue
# tissues_df.loc[tissues_df.tissue == 'thyroid_gland', 'tissue'] = 'thyroid'

# Set all subtissues that don't have both tumor and normal to 'other_tissues'
# tissues_df.loc[~tissues_df.tissue.isin(tissues), 'tissue'] = 'other_tissues'

# Print all tissues with corresponding number of tumor and normal subtissues
#print('Tissues:', tissues_df.groupby('tissue').subtissue.nunique())

#print(tissues_df.head())
#print(tissues_df['path'].head())

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
    outdir = os.path.join('/home/lnemati/pathway_crosstalk/results/breakpoints', filename.replace('.csv', ''), metric)
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

print("Done: breakpoints.py")


# OLD CODE: no multiprocessing
#
#for metric in metrics:
#    for i, row in tissues_df.iterrows():
#        tissue = row.tissue
#        path = row.path
#        subtissue = row.subtissue
#        condition = row.condition
#        path = row.path
#
#        print('tissue:', tissue)
#        print('subtissue:', subtissue)
#
#        # Make a copy of all_interactions for this tissue
#        tissue_interactions = all_interactions.copy()
#        tissue_interactions[metric] = 0
#
#        # Read the interactions table for this subtissue
#        df = pd.read_csv(os.path.join(path, 'interactions', filename), index_col='interaction')
#        # Fill nan values with 0
#        df = df.fillna(0)
#
#        # Add the metric values to the tissue_interactions dataframe
#        tissue_interactions[metric] = df[metric]
#
#        # Add avg values to normal_vals and tumor_vals, use nans for missing
#        # Rows are interactions, columns are tissues
#        if condition == 'normal':
#            normal_vals[metric][subtissue] = np.nan
#            normal_vals[metric].loc[tissue_interactions.index, subtissue] = tissue_interactions[metric]
#        elif condition == 'tumor':
#            tumor_vals[metric][subtissue] = np.nan
#            tumor_vals[metric].loc[tissue_interactions.index, subtissue] = tissue_interactions[metric]
#
#        # For each tissue count how many interactions have missing genes and add it to the dataframe
#        missing_genes[subtissue] = df['missing_genes']
#
#    # Make output dir with the name of the file
#    outdir = os.path.join('/home/lnemati/pathway_crosstalk/results/breakpoints', filename.replace('.csv', ''), metric)
#    os.makedirs(outdir, exist_ok=True)
#
#    # Save the values per tissue to a csv file
#    print('Saving: ', f'{outdir}/normal_values.csv')
#    normal_vals[metric].to_csv(f'{outdir}/normal_values.csv')
#    print('Saving: ', f'{outdir}/tumor_values.csv')
#    tumor_vals[metric].to_csv(f'{outdir}/tumor_values.csv')
#    print('Saving: ', f'{outdir}/missing_genes.csv')
#    missing_genes.to_csv(f'{outdir}/missing_genes.csv')
