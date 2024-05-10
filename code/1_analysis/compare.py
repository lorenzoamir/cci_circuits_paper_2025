import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare results of individual WGCNA analysis')

parser.add_argument('--dir_list', type=str, help='List of directories to compare, delimited by commas')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results'

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]

check_files = [
    'general_info.txt',
    'interactions_info.txt',
    'degree_df.csv.gz',
    'stats.txt',
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

df = pd.DataFrame(
    columns=[
        'name',
        'condition',
        'n_samples',
        'n_genes',
        'n_modules',
        'avg_module_size',
        'avg_degree',
        'avg_intramodular_degree',
        'tot_interactions',
        'n_interactions_same_module',
        'interactions_auroc',
        'interactions_ranksum_U',
        'interactions_ranksum_p',
    ]
)

def get_name_condition(d):
    if '/tumor/' in d: 
        condition = 'tumor'
    elif '/normal/' in d:
        condition = 'normal'
    else:
        print('Error: Directory {} does not contain /tumor or /normal'.format(d), file=sys.stderr)
        print('Error: Directory {} does not contain /tumor or /normal'.format(d), file=sys.stdout)
        return None, None
    
    name = d.split('/')[-1]
    return name, condition

for d in dir_list:

    print(d)
    name, condition = get_name_condition(d)
    print(name)
    print(condition)
    
    general_info = os.path.join(d, 'general_info.txt')
    interactions_info = os.path.join(d, 'interactions_info.txt')
    degree_df = os.path.join(d, 'degree_df.csv.gz')
    stats = os.path.join(d, 'stats.txt')

    # add name to df index
    df.loc[name, 'name'] = name
    df.loc[name, 'condition'] = condition

    # Read general info
    print('Reading general info')
    with open(general_info, 'r') as f:
        for line in f:
            key = line.strip().split(': ')[0]
            value = line.strip().split(': ')[1:]
            if len(value) == 1:
                value = value[0]
            if key == 'n_samples':
                df.loc[name, 'n_samples'] = int(value)
            elif key == 'n_genes':
                df.loc[name, 'n_genes'] = int(value)
            elif key == 'n_modules':
                df.loc[name, 'n_modules'] = int(value)
            elif key == 'module_sizes':
                # join back values with ', ' as separator
                sizes = ', '.join(value)
                # remove { and } from string
                sizes = sizes.replace('{', '')
                sizes = sizes.replace('}', '')
                # make list by separating by ', '
                sizes = sizes.split(', ')
                sizes = sizes[1::2]
                sizes = [int(x) for x in sizes]
                df.loc[name, 'avg_module_size'] = np.mean(sizes)
            
    # Read degree_df
    print('Reading degree_df')
    degree_df = pd.read_csv(degree_df, index_col=0)
    avg_degree = degree_df['degree'].mean()
    avg_intramodular_degree = degree_df['intramodular_degree'].mean()
    df.loc[name, 'avg_degree'] = avg_degree
    df.loc[name, 'avg_intramodular_degree'] = avg_intramodular_degree

    # Read interactions_info
    print('Reading interactions_info')
    with open(interactions_info, 'r') as f:
        for line in f:
            key = line.strip().split(': ')[0]
            value = line.strip().split(': ')[1:]
            if len(value) == 1:
                value = value[0]
            if key == 'total_interactions':
                df.loc[name, 'tot_interactions'] = int(value)
            elif key == 'n_interactions_same_module':
                df.loc[name, 'n_interactions_same_module'] = int(value)

    # Read stats
    print('Reading stats')
    with open(stats, 'r') as f:
        for line in f:
            key = line.strip().split(': ')[0]
            value = line.strip().split(': ')[1:]
            if len(value) == 1:
                value = value[0]
            if key == 'auroc':
                df.loc[name, 'interactions_auroc'] = float(value)
            elif key == 'ranksums_U_all':
                df.loc[name, 'interactions_ranksum_U'] = float(value)
            elif key == 'ranksums_p_all':
                df.loc[name, 'interactions_ranksum_p'] = float(value)
    
# Save results
df.to_csv(os.path.join(output, 'compare_results.csv'), index=False)

print('Done: compare.py')

