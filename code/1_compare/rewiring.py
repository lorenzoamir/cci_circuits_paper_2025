import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find LR pairs that are consistently rewired')

parser.add_argument('-r', '--lr_resource', type=str, help='Path to LR resource file')
parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

# Set defaults
if args.lr_resource is None:
    args.lr_resource = '/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/consensus.csv'
if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

# Load the LR resource file
print("Loading LR resource file: ", args.lr_resource)
lr_resource = pd.read_csv(args.lr_resource)
lr_resource = lr_resource[['ligand', 'receptor']]
print(lr_resource.head())
print()

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

t_same_df = lr_resource.copy()
n_same_df = lr_resource.copy()

def rewiring(path, outdf, lr_resource):
    # Get cancer type from directory name
    cancer_type = os.path.basename(os.path.dirname(path))

    df = pd.read_csv(path)
    outdf[cancer_type] = False

    # Set to True all ligand-receptor pairs that are in df
    outdf.loc[outdf['ligand'].isin(df['ligand']) & outdf['receptor'].isin(df['receptor']), cancer_type] = True

    return outdf
    

for t_same in t_same_n_diff:
    t_same_df = rewiring(t_same, t_same_df, lr_resource)

# Count tissues for each LR pair
t_same_df['total'] = t_same_df.iloc[:, 2:].sum(axis=1)

# columns order: ligand, receptor, total, rest of the columns
cols = t_same_df.columns.tolist()
cols = cols[:2] + cols[-1:] + cols[2:-1]

t_same_df = t_same_df[cols]

# Sort by total
t_same_df = t_same_df.sort_values(by='total', ascending=False)

# Save the t_same_df
t_same_df.to_csv(args.outputdir + 't_same_n_diff.csv', index=False)

for n_same in n_same_t_diff:
    n_same_df = rewiring(n_same, n_same_df, lr_resource)

# Count tissues for each LR pair
n_same_df['total'] = n_same_df.iloc[:, 2:].sum(axis=1)

# columns order: ligand, receptor, total, rest of the columns
cols = n_same_df.columns.tolist()
cols = cols[:2] + cols[-1:] + cols[2:-1]

n_same_df = n_same_df[cols]

# Sort by total
n_same_df = n_same_df.sort_values(by='total', ascending=False)

# Save the n_same_df
n_same_df.to_csv(args.outputdir + 'n_same_t_diff.csv', index=False)

print('Done: rewiring.py')
