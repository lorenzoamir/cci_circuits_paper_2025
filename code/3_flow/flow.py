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

# Print each path and whether it exists
for path in tissues_df['path']:
    print(path, os.path.exists(path))

print('Tissues dataframe:', tissues_df.shape)
print(tissues_df.head())

# Find tissues that have both tumor and normal conditions
tissues = tissues_df.groupby('tissue').condition.nunique()
tissues = tissues[tissues == 2].index
print('Tissues with both tumor and normal:', len(tissues))

# Remove normal testis subtissues
tissues_df = tissues_df[~((tissues_df.tissue == 'testis') & (tissues_df.condition == 'normal'))]
tissues = [x for x in tissues if x != 'testis']

# Join thyroid_gland and thyroid subtissues into the same tissue
tissues_df.loc[tissues_df.tissue == 'thyroid_gland', 'tissue'] = 'thyroid'

# Set all subtissues that don't have both tumor and normal to 'other_tissues'
tissues_df.loc[~tissues_df.tissue.isin(tissues), 'tissue'] = 'other_tissues'

# Print all tissues with corresponding number of tumor and normal subtissues
print('Tissues:', tissues_df.groupby('tissue').subtissue.nunique())

print(tissues_df.head())
print(tissues_df['path'].head())

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

interactions_files = ['ccc_lr_pairs.csv', 'intact_direct.csv', 'intact_physical.csv', 'intact_association.csv']
metrics = ['same_module', 'adj', 'corr', 'kme_corr'] #, 'min_adj', 'mean_adj', 'min_kme_corr', 'mean_kme_corr']

for filename in interactions_files:
    print('filename:', filename)

    all_interactions = pd.read_csv(os.path.join('/home/lnemati/pathway_crosstalk/data/interactions', filename), index_col='interaction')
    
    # Init columns for each metric
    for metric in metrics:
        all_interactions['avg_normal'] = 0
        all_interactions['avg_tumor'] = 0

    tissues_cat = {k: [] for k in metrics}
    normal_cat = {k: [] for k in metrics}
    tumor_cat = {k: [] for k in metrics}
    outcome_cat = {k: [] for k in metrics}
    counts = {k: [] for k in metrics}
    colors = {k: [] for k in metrics}

    normal_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}
    tumor_vals = {k: pd.DataFrame(index=all_interactions.index) for k in metrics}

    for metric in metrics:
        for tissue in tissues_df.tissue.unique():
            tissue_name = tissue
            print('tissue:', tissue_name)
       
            tissue_interactions = all_interactions.copy() 
            tissue_interactions[metric] = 0

            # Find all interactions.csv files in the normal and tumor directories
            normal_dfs = [os.path.join(n_dir, 'interactions', filename) for n_dir in normal_dirs[tissue_name]]   
            tumor_dfs = [os.path.join(t_dir, 'interactions', filename) for t_dir in tumor_dirs[tissue_name]]

            # Read the dataframes
            normal_dfs = [pd.read_csv(df, index_col='interaction') for df in normal_dfs]
            tumor_dfs = [pd.read_csv(df, index_col='interaction') for df in tumor_dfs]

            # Print shapes
            print('Normal shapes:', [df.shape for df in normal_dfs])
            print('Tumor shapes:', [df.shape for df in tumor_dfs])

            # Merge normal and tumor interactions with the full interaction network
            print('Merging')
            new_dfs = []
            for df in normal_dfs:
                n_tmp = tissue_interactions.copy()
                n_tmp['metric'] = 0
                n_tmp.loc[df.index, metric] = df[metric] # Add the values of the metric to the interaction network
                new_dfs.append(n_tmp) # Add the new dataframe to the list
            normal_dfs = new_dfs

            new_dfs = []
            for df in tumor_dfs:
                t_tmp = tissue_interactions.copy()
                t_tmp['metric'] = 0
                t_tmp.loc[df.index, metric] = df[metric] # Add the values of the metric to the interaction network
                new_dfs.append(t_tmp) # Add the new dataframe to the list
            tumor_dfs = new_dfs
           
            # Make sure all dfs follow the same order
            normal_dfs = [df.loc[tissue_interactions.index] for df in normal_dfs]
            tumor_dfs = [df.loc[tissue_interactions.index] for df in tumor_dfs]

            print('Normal shapes:', [df.shape for df in normal_dfs])
            print('Tumor shapes:', [df.shape for df in tumor_dfs])

            # Average the value of the metric in the normal and tumor networks on all subtissues
            tissue_interactions['avg_normal'] = sum([df[metric] for df in normal_dfs]) / len(normal_dfs)
            tissue_interactions['avg_tumor'] = sum([df[metric] for df in tumor_dfs]) / len(tumor_dfs)

            # Add avg values to normal_vals and tumor_vals, use nans for missing
            # Rows are interactions, columns are tissues
            #normal_vals[metric][tissue_name] = np.nan
            #tumor_vals[metric][tissue_name] = np.nan
            normal_vals[metric].loc[tissue_interactions.index, tissue_name] = tissue_interactions['avg_normal']
            tumor_vals[metric].loc[tissue_interactions.index, tissue_name] = tissue_interactions['avg_tumor']

            # Get flow for each outcome
            tissues_cat[metric] += [format_string(tissue_name, newline=False)] * 3
            normal_cat[metric] += ['Same<br>Module', 'Same<br>Module', 'Different<br>Modules']
            tumor_cat[metric] += ['Same<br>Module', 'Different<br>Modules', 'Same<br>Module']
            outcome_cat[metric] += ['Both', 'Normal<br>Only', 'Tumor<br>Only']
            
            # Both is the fraction of tissues in which the interaction is in the same module in both normal and tumor
            # So its the minimum of the two fractions
            both = tissue_interactions[['avg_normal', 'avg_tumor']].min(axis=1)
            normal_only = tissue_interactions['avg_normal'] - both
            tumor_only = tissue_interactions[f'avg_tumor'] - both

            # Sum values
            both = both.sum()
            normal_only = normal_only.sum()
            tumor_only = tumor_only.sum()

            counts[metric] += [both, normal_only, tumor_only]
            colors[metric] += [graycolor2, ncolor, tcolor]
            
            # Save fractions to all_interactions
            #all_interactions.loc[tissue_interactions.index, f'both_score'] += both
            all_interactions.loc[tissue_interactions.index, 'avg_normal'] += tissue_interactions['avg_normal'] / tissues_df.tissue.nunique()
            all_interactions.loc[tissue_interactions.index, 'avg_tumor'] += tissue_interactions['avg_tumor'] / tissues_df.tissue.nunique()

        # After reading all tissues and make the actual plots
        # Normalize counts
        counts[metric] = np.array(counts[metric])

        # Save categories and counts as a csv file
        df = pd.DataFrame({'Tissue': tissues_cat[metric], 'Normal': normal_cat[metric], 'Tumor': tumor_cat[metric], 'Outcome': outcome_cat[metric], 'Counts': counts[metric]})
        # Make output dir with the name of the file
        outdir = os.path.join('/home/lnemati/pathway_crosstalk/results/flow', filename.replace('.csv', ''), metric)
        os.makedirs(outdir, exist_ok=True)
        print('Saving: ', f'{outdir}/flow_diagram_data.csv')
        df.to_csv(f'{outdir}/flow_diagram_data.csv', index=False)
       
        # Save the values per tissue to a csv file
        print('Saving: ', f'{outdir}/normal_values.csv')
        normal_vals[metric].to_csv(f'{outdir}/normal_values.csv')
        print('Saving: ', f'{outdir}/tumor_values.csv')
        tumor_vals[metric].to_csv(f'{outdir}/tumor_values.csv')

        # Remove interactions that only have 0 values
        # keep = all_interactions[(all_interactions[f'both_score'] > 0) | (all_interactions[f'normal_score'] > 0) | (all_interactions[f'tumor_score'] > 0)].index
        # all_interactions = all_interactions.loc[keep]
        all_interactions[f'diff'] = all_interactions['avg_tumor'] - all_interactions['avg_normal']

        # Sort by difference between tumor and normal
        sorted_interactions = all_interactions.sort_values('diff', ascending=False)
        sorted_interactions = sorted_interactions[['avg_tumor', 'avg_normal', 'diff']]
        print('Saving: ', f'{outdir}/interactions_with_scores.csv')
        sorted_interactions.to_csv(f'{outdir}/interactions_with_scores.csv')
        
        if metric == 'same_module':

            # Make parallel categories plot, use colors to distinguish between normal and tumor
            fig = go.Figure(data=[go.Parcats(
                dimensions=[
                    {'label': 'Tissue',
                     'values': tissues_cat[metric]},
                    {'label': 'Normal',
                     'values': normal_cat[metric]},
                    {'label': 'Tumor',
                     'values': tumor_cat[metric]},
                    {'label': 'Outcome',
                     'values': outcome_cat[metric]}
                ],
                counts=counts[metric],
                line={'shape': 'hspline', 'color': colors[metric]},
                )]
            )

            # Larger margins, larger font, black font color
            fig.update_layout(
                #title_text='Flow of LR pairs',
                margin=dict(l=250, r=120, t=50, b=50),
                font_size=45,
                font_color='black'
            )

            # Save figure make the figure big
            print('Saving image to:', outdir)
            fig.write_image(os.path.join(outdir, "flow.pdf"), width=3000, height=1800)
            fig.write_image(os.path.join(outdir, "flow.png"), width=3000, height=1800)

            print()

print("Done: flow.py")
