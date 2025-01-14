import pandas as pd
import numpy as np
import gseapy as gp
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Get average correlation for ccc and intact')

parser.add_argument('--dir-list', type=str, help='List of directories to compare, delimited by commas')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results/comparison'
# Create output directory if it does not exist
os.makedirs(output, exist_ok=True)

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]
# Add interactions to the paths
dir_list = [os.path.join(d, 'interactions') for d in dir_list]

check_files = [
    'ccc_lr_pairs.csv',
    'intact_direct.csv',
    'intact_physical.csv',
    'intact_association.csv',
]

# Sanity check
for d in dir_list:
    if not os.path.exists(d):
        # error message to standard error and standard output
        print('Error: Directory {} does not exist'.format(d), file=sys.stderr)
        print('Error: Directory {} does not exist'.format(d), file=sys.stdout)
    # Check if all files are present
    for f in check_files:
        if not os.path.exists(os.path.join(d, f)):
            print('Error: File {} does not exist in directory {}'.format(f, d), file=sys.stderr)
            print('Error: File {} does not exist in directory {}'.format(f, d), file=sys.stdout)

def get_name_condition(d):
    subtissue = d.split('/')[-2]
    condition = d.split('/')[-3]
    major_tissue = d.split('/')[-4]
    
    return subtissue, condition, major_tissue

all_genes = set()
all_interactors = set()

normal_tissues = []
tumor_tissues = []

# Init dataframe, index will be tissues, columns will be average corr and adj and median corr and adj
df = pd.DataFrame(columns=['subtissue', 'condition', 'major_tissue', 'interaction_type', 'value_type', 'value'])

print('Reading all files')
for d in dir_list:

    print(d)
    subtissue, condition, major_tissue = get_name_condition(d)
    print(subtissue, condition, major_tissue)
   
    for file in check_files:
        int_df = pd.read_csv(os.path.join(d, file))
        # Fill nans in corr and adj with 0

        rows = [
            pd.DataFrame({ 
                'subtissue': subtissue,
                'condition': condition,
                'major_tissue': major_tissue,
                'interaction_type': file ,
                'value_type': 'average_corr',
                'value': int_df['corr'].mean()
            }, index=[0]),
            pd.DataFrame({
                'subtissue': subtissue,
                'condition': condition,
                'major_tissue': major_tissue,
                'interaction_type': file,
                'value_type': 'average_adj',
                'value': int_df['adj'].mean()
            }, index=[0]),
            pd.DataFrame({
                'subtissue': subtissue,
                'condition': condition,
                'major_tissue': major_tissue,
                'interaction_type': file,
                'value_type': 'median_corr',
                'value': int_df['corr'].median()
            }, index=[0]),
            pd.DataFrame({
                'subtissue': subtissue,
                'condition': condition,
                'major_tissue': major_tissue,
                'interaction_type': file,
                'value_type': 'median_adj',
                'value': int_df['adj'].median()
            }, index=[0])
        ]

        df = pd.concat([df, *rows], ignore_index=True)

df.to_csv(os.path.join(output, 'averages.csv'), index=False)

print('Done: averages.py')
