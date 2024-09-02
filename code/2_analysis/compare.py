import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare results of individual WGCNA analysis')

parser.add_argument('--dir-list', type=str, help='List of directories to compare, delimited by commas')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results'

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]

check_files = [
    'general_info.txt',
    'interactions_info.txt',
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
        'n_modules',
        'n_genes',
        'n_samples',
        'sft_power',
        'avg_module_size',
        'avg_coev_tom_correlation',
        'tot_interactions',
        'n_interactions_same_module',
        'interactions_auroc',
        'interactions_mannwhitneyu_U',
        'interactions_mannwhitneyu_p',
        'avg_degree',
        'avg_intramodular_degree',
        'var_degree',
        'var_intramodular_degree',
        'powerlaw_alpha',
        'powerlaw_xmin',
        'powerlaw_sigma',
        'avg_path_length',
        'diameter',
        'clustering_global',
        'modularity',
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

def remove_duplicates(file):
    # Each line of the file follows "key: value"
    # If a key appears more than once, only the first value is kept
    with open(file, 'r') as f:
        lines = f.readlines()
    with open(file, 'w') as f:
        seen = set()
        for line in lines:
            if line.split(': ')[0] not in seen:
                f.write(line)
                seen.add(line.split(': ')[0])

for d in dir_list:

    print(d)
    name, condition = get_name_condition(d)
    print(name)
    print(condition)
    
    general_info = os.path.join(d, 'general_info.txt')
    
    # Remove duplicates from general_info
    #remove_duplicates(general_info)

    interactions_info = os.path.join(d, 'interactions_info.txt')
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
            elif key == 'sft_power':
                df.loc[name, 'sft_power'] = float(value)
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
            elif key == 'average_tom_coev_correlation':
                df.loc[name, 'avg_coev_tom_correlation'] = float(value)
            elif key == 'avg_degree':
                df.loc[name, 'avg_degree'] = float(value)
            elif key == 'avg_intramodular_degree':
                df.loc[name, 'avg_intramodular_degree'] = float(value)
            elif key == 'var_degree':
                df.loc[name, 'var_degree'] = float(value)
            elif key == 'var_intramodular_degree':
                df.loc[name, 'var_intramodular_degree'] = float(value)
            elif key == 'powerlaw_alpha':
                df.loc[name, 'powerlaw_alpha'] = float(value)
            elif key == 'powerlaw_xmin':
                df.loc[name, 'powerlaw_xmin'] = float(value)
            elif key == 'powerlaw_sigma':
                df.loc[name, 'powerlaw_sigma'] = float(value)
            elif key == 'avg_path_length':
                df.loc[name, 'avg_path_length'] = float(value)
            elif key == 'diameter':
                df.loc[name, 'diameter'] = float(value)
            elif key == 'clustering_global':
                df.loc[name, 'clustering_global'] = float(value)
            elif key == 'modularity':
                df.loc[name, 'modularity'] = float(value)
                
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
            elif key == 'mannwhitneyu_U_all':
                df.loc[name, 'interactions_mannwhitneyu_U'] = float(value)
            elif key == 'mannwhitneyu_p_all':
                df.loc[name, 'interactions_mannwhitneyu_p'] = float(value)
    
# Save results
df.to_csv(os.path.join(output, 'compare_results.csv'), index=False)

print('Done: compare.py')

