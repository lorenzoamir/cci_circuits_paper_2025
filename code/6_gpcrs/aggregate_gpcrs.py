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

data_dir='/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal/'

# Save tissues_df to a csv file
tissues_df = pd.read_csv('/home/lnemati/pathway_crosstalk/data/tissues.csv')

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

df = pd.DataFrame(columns=['gpcr', 'interactor', 'condition','interactor_type', 'tissue', 'subtissue', 'value_type', 'value'])

# For each tissue read all corresponding files and average the values
for i, row in tissues_df.iterrows():
    tissue = row.tissue
    subtissue = row.subtissue
    path = row.path

    gpcr_df = pd.read_csv(os.path.join(path, 'gpcrs.csv'))

    # Add tissue and subtissue columns
    gpcr_df['tissue'] = tissue
    gpcr_df['subtissue'] = subtissue

    # Concatenate to the main dataframe
    df = pd.concat([df, gpcr_df])

# Save the dataframe to a csv file
df.to_csv('/home/lnemati/pathway_crosstalk/results/gpcrs/all_gpcr_interactions.csv', index=False)

print('Done: aggregate_gpcrs.py')
