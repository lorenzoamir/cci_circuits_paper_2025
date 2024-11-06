import pandas as pd
import os
import numpy as np
import PyWGCNA
import gseapy as gp
import argparse
import plotly.graph_objects as go
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

def format_string(string, newline=False):
    '''
    replaces _ with space and capitalizes the first letter of each word.
    finds which of the spaces can be replaced with a newline as to make the
    two resulting lines of text as close to equal length as possible.
    '''
    # Replace _ with space
    string = string.replace('_', ' ').title()
    
    if newline:
        # If string has less than 20 characters, don't bother
        if len(string) < 10:
            return string
        # Find the space that makes the two lines of text as close to equal length as possible
        diffs = []

        for i, words in enumerate(string.split()):
            diffs.append(abs(len(' '.join(string.split()[:i])) - len(' '.join(string.split()[i:]))))
            
        # Find which space corresponds to the minimum difference and replace it with <br>
        n_space = diffs.index(min(diffs)) # e.g. the second ' ' in the string
        string = ' '.join(string.split()[:n_space]) + '<br>' + ' '.join(string.split()[n_space:])

    return string

filename = 'all_ccc_gene_pairs.csv'
metrics = ['adj', 'kme_corr', 'corr']

print('filename:', filename)

all_interactions = pd.read_csv(os.path.join('/home/lnemati/pathway_crosstalk/data/interactions', filename), index_col='interaction')

# Init columns for each metric
for metric in metrics:
    all_interactions[f'{metric}_normal'] = 0
    all_interactions[f'{metric}_tumor'] = 0

normal_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}
tumor_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}
missing_genes = pd.DataFrame(index=all_interactions.index)

for metric in metrics:
    for i, row in tissues_df.iterrows():
        tissue = row.tissue
        path = row.path
        subtissue = row.subtissue
        condition = row.condition
        path = row.path

        print('tissue:', tissue)
        print('subtissue:', subtissue)

        # Make a copy of all_interactions for this tissue
        tissue_interactions = all_interactions.copy() 
        tissue_interactions[metric] = 0
    
        # Read the interactions table for this subtissue
        df = pd.read_csv(os.path.join(path, 'interactions', filename), index_col='interaction')
        # Fill nan values with 0
        df = df.fillna(0)
        
        # Add the metric values to the tissue_interactions dataframe
        tissue_interactions[metric] = df[metric]
        
        # Add avg values to normal_vals and tumor_vals, use nans for missing
        # Rows are interactions, columns are tissues
        if condition == 'normal':
            normal_vals[metric][subtissue] = np.nan
            normal_vals[metric].loc[tissue_interactions.index, subtissue] = tissue_interactions[metric]
        elif condition == 'tumor':
            tumor_vals[metric][subtissue] = np.nan
            tumor_vals[metric].loc[tissue_interactions.index, subtissue] = tissue_interactions[metric]
        else:
            print('Error: condition not tumor or normal')
            print('Error: condition not tumor or normal', file=sys.stderr)
            sys.exit(1)

        # For each tissue count how many interactions have missing genes and add it to the dataframe
        missing_genes[subtissue] = df['missing_genes']

    # Make output dir with the name of the file
    outdir = os.path.join('/home/lnemati/pathway_crosstalk/results/breakpoints', filename.replace('.csv', ''), metric)
    os.makedirs(outdir, exist_ok=True)
   
    # Save the values per tissue to a csv file
    print('Saving: ', f'{outdir}/normal_values.csv')
    normal_vals[metric].to_csv(f'{outdir}/normal_values.csv')
    print('Saving: ', f'{outdir}/tumor_values.csv')
    tumor_vals[metric].to_csv(f'{outdir}/tumor_values.csv')
    print('Saving: ', f'{outdir}/missing_genes.csv')
    missing_genes.to_csv(f'{outdir}/missing_genes.csv')
        
print("Done: breakpoints.py")
