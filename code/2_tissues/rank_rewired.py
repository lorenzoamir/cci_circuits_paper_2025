import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find LR pairs that are consistently rewired')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

lr_resource = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')

# Searching for t_same_n_diff and n_same_t_diff files in subdirectories
parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'
t_same_n_diff = []
n_same_t_diff = []
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.endswith('t_same_n_diff.csv'):
            t_same_n_diff.append(os.path.join(root, file))
        elif file.endswith('n_same_t_diff.csv'):
            n_same_t_diff.append(os.path.join(root, file))
print('Found the following t_same_n_diff files: ')
print(t_same_n_diff)
print()
print('Found the following n_same_t_diff files: ')
print(n_same_t_diff)
print()

print('lr resource: ')
print(lr_resource.head())

t_same_df = lr_resource[['interaction']].copy()
n_same_df = lr_resource[['interaction']].copy()

print('t_same_df: ')
print(t_same_df.head())

print('n_same_df: ')
print(n_same_df.head())

def rewiring(path, outdf, lr_resource):
    # Get cancer type from directory name
    cancer_type = os.path.basename(os.path.dirname(path))

    df = pd.read_csv(path)

    outdf[cancer_type] = False

    # Set to True all interactions that are in df
    outdf.loc[outdf['interaction'].isin(df['interaction']), cancer_type] = True
    
    return outdf

# ------ T SAME N DIFF ------

for t_same in t_same_n_diff:
    print(t_same)
    try:
        t_same_df = rewiring(t_same, t_same_df, lr_resource)
    except:
        print('Error: ', t_same)
        # Remove cancer type col from t_same_df
        cancer_type = os.path.basename(os.path.dirname(t_same))
        t_same_df = t_same_df.drop(columns=[cancer_type])

print()
print('t_same_df: ')
print(t_same_df.head())

# Count tissues for each LR pair
t_same_df['total'] = t_same_df.iloc[:, 1:].sum(axis=1)

# columns order: interaction, total, rest of the columns 
cols = t_same_df.columns.tolist()
cols = ['interaction', 'total'] + [col for col in cols if col not in ['interaction', 'total']]

t_same_df = t_same_df[cols]

# Sort by total
t_same_df = t_same_df.sort_values(by='total', ascending=False)

# Save the t_same_df
t_same_df.to_csv(args.outputdir + 't_same_n_diff.csv', index=False)

# ------ N SAME T DIFF ------

for n_same in n_same_t_diff:
    print(n_same)
    try:
        n_same_df = rewiring(n_same, n_same_df, lr_resource)
    except:
        print('Error: ', n_same)
        # Remove cancer type col from n_same_df
        cancer_type = os.path.basename(os.path.dirname(n_same))
        n_same_df = n_same_df.drop(columns=[cancer_type])

# Count tissues for each LR pair
n_same_df['total'] = n_same_df.iloc[:, 1:].sum(axis=1)

# columns order: interaction, total, rest of the columns
cols = n_same_df.columns.tolist()
cols = ['interaction', 'total'] + [col for col in cols if col not in ['interaction', 'total']]
n_same_df = n_same_df[cols]

# Sort by total
n_same_df = n_same_df.sort_values(by='total', ascending=False)

# Save the n_same_df
n_same_df.to_csv(args.outputdir + 'n_same_t_diff.csv', index=False)

print('Done: rank_rewired.py')
