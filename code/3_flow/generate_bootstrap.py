import pandas as pd
import os
import numpy as np
import PyWGCNA
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description='Generate bootstrap samples')

# list of tumor and normal directories
parser.add_argument('--tissues', type=str, help='List of tissue directories separated by space', required=True)

args = parser.parse_args()

output_dir = '/home/lnemati/pathway_crosstalk/results/flow/bootstrap'

# Inside each tissue directory, there are two directories: normal and tumor. Find all subdirectories of normal and tumor
print('Input tissues:', args.tissues)

tumor_dirs = {}
normal_dirs = {}

all_tissues = args.tissues.split(" ")

# Remove empty strings and trailing characters
all_tissues = [x.strip() for x in all_tissues if x.strip() != '']

for tissue_dir in all_tissues:
    # Get tissue name from path
    tissue_name = tissue_dir.split("/")[-1]

    # Get tumor subdirectories
    tumor_dirs[tissue_name] = [os.path.join(tissue_dir, 'tumor', x) for x in os.listdir(os.path.join(tissue_dir, 'tumor'))]
    # Get normal subdirectories
    normal_dirs[tissue_name] = [os.path.join(tissue_dir, 'normal', x) for x in os.listdir(os.path.join(tissue_dir, 'normal'))]

tissues_cat = []
normal_cat = []
tumor_cat = []
outcome_cat = []
counts = []
colors = []

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

# Read the full interaction network
print('Reading interaction network')
all_interactions = pd.read_csv('/projects/bioinformatics/DB/CellCellCommunication/WithEnzymes/cpdb_cellchat_enz.csv')

# Extract genes for complex A and B for each row
complex_a_genes = all_interactions.apply(lambda row: [row[f'interactor{i}'] for i in range(1, row['num_interactors_a'] + 1) if pd.notna(row[f'interactor{i}'])], axis=1)
complex_b_genes = all_interactions.apply(lambda row: [row[f'interactor{i}'] for i in range(row['num_interactors_a'] + 1, row['num_interactors_a'] + row['num_interactors_b'] +     1) if pd.notna(row[f'interactor{i}'])], axis=1)

complex_a_genes.index = all_interactions['interaction']
complex_b_genes.index = all_interactions['interaction']

complex_a = list(complex_a_genes.apply(lambda x: '+'.join([gene for gene in x]) + '_'))
complex_b = list(complex_b_genes.apply(lambda x: '+'.join([gene for gene in x])))

# Generate all possible combinations using MultiIndex and concatenate strings
print('Generating all possible combinations')
bootstrap = pd.MultiIndex.from_product([all_interactions['interaction'], all_interactions['interaction']])

# Convert to DataFrame for a cleaner format (optional)
print('Converting to DataFrame')
bootstrap = pd.DataFrame(
    complex_a_genes.loc[bootstrap.get_level_values(0)].values + complex_b_genes.loc[bootstrap.get_level_values(1)].values,
    columns=['all_genes']
)
bootstrap['n_genes'] = bootstrap['all_genes'].apply(lambda x: len(x))

N = 5000

df_all = pd.DataFrame()
# Bootstrap
print('Sampling')
for n_genes in range(2,8):
    print(n_genes)
    df = bootstrap.query('n_genes == @n_genes')
    random = np.random.choice(range(len(df)), N, replace=True)
    df = df.iloc[random]
    
    # Add to bootstrap_all
    df_all = pd.concat([df_all, df])

    # Remove n_genes column
    df.drop('n_genes', axis=1, inplace=True)

    # Save dataframe of bootstrap samples
    df.to_csv(os.path.join(output_dir, f'bootstrap_{n_genes}.csv'), index=False)

# Save one bootstrap_all dataframe containing merged bootstrap datasets
df_all.to_csv(os.path.join(output_dir, 'bootstrap_all.csv'), index=False)

print("Done: generata_bootstrap.py")
