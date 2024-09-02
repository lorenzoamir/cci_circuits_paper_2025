import pandas as pd
import numpy as np
import re
import os
import multiprocessing
from tqdm import tqdm
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
from itertools import combinations

def extract_gene_symbols(mutation_string):
    if pd.isna(mutation_string) or mutation_string == '':
        return []

    # Split the string into individual entries
    gene_entries = mutation_string.split(',')

    # Filter entries to only those with parentheses
    allowed_strings = ['(p.', '(amp', '(del']

    filtered_entries = [entry for entry in gene_entries if any([ x in entry for x in allowed_strings]) and ')' in entry]

    # Remove the parentheses and any content within them
    cleaned_entries = [re.sub(r'\([^()]*\)', '', entry).strip() for entry in filtered_entries]

    # Return the cleaned entries
    return cleaned_entries

#gene_names = pd.read_csv('/home/lnemati/resources/biomart/ensembl_to_symbol_filtered.csv.gz')['Gene name'].values
output_dir = '/home/lnemati/pathway_crosstalk/results/coevolution/tumor'
df_complete = pd.read_csv('/projects/bioinformatics/DB/CancerTracer/intra_heterogeneity_data.txt', sep='\t')
df_complete = df_complete.set_index('Sample', drop=False)
df_complete = df_complete.replace({'-': ''})

# Remove wrongly formatted row
df_complete = df_complete[~df_complete.index.isin(['p075_p3'])]

df_complete['Trunk'] = df_complete['Trunk_mutation'].apply(extract_gene_symbols)
df_complete['Branch'] = df_complete['Branch_mutation'].apply(extract_gene_symbols)
df_complete['Private'] = df_complete['Private_mutation'].apply(extract_gene_symbols)

print('Subsetting to only coding, amp and del')
keep = df_complete.apply(lambda x: any([len(x['Trunk']), len(x['Branch']), len(x['Private'])]) != 0, axis=1).values
df_complete = df_complete[keep]
print('Subsetted shape:', df_complete.shape)

mutation_types = ['Trunk', 'Branch', 'Private', 'Patient']

# Initialize the modules
all_patients = set(df_complete.index)

# Change Cancertype to lower case, replace spaces with underscores
df_complete['Cancertype'] = df_complete['Cancertype'].str.lower().str.replace(' ', '_')
# Add cancertype to patient names
df_complete.index = df_complete['Cancertype'] + '_' + df_complete.index

# Duplicate dataframe
df_copy = df_complete.copy()
# Set cancer type to 'All Cancers'
df_copy['Cancertype'] = 'all_cancers'

# Concatenate the two dataframes
print('Original shape:', df_complete.shape)
print('Copy shape:', df_copy.shape)
df_complete = pd.concat([df_complete, df_copy])
print('Concatenated shape:', df_complete.shape)

print('Running for each cancertype')
num_processes = multiprocessing.cpu_count()

def process_cancer_type(df):
    modules = {}

    if df['Cancertype'].nunique() > 1:
        # Raise and error if there are multiple cancertypes
        raise ValueError('Multiple cancertypes in dataframe: ', df['Cancertype'].unique())
    
    cancer_type = df['Cancertype'].iloc[0]

    group = df.copy()
    group = group[['Trunk', 'Branch', 'Private']]
    group = group.stack().reset_index()

    group = group.rename(columns={'level_0': 'patient', 'level_1': 'mutation', 0: 'genes'})
    #group = group.rename(columns={'level_1': 'mutation', 0: 'genes', 'Sample': 'patient'})

    for col in group.columns:
        print(col)
        print(group.head(3)[col])

    for col in group.columns:
        print(col, group[col].dtype)

    group['module'] = group['patient'] + '_' + group['mutation']
    group = group[['genes', 'module', 'patient', 'mutation']]

    for mutation_type in mutation_types:
        print(mutation_type)
        df = group.copy()

        if mutation_type in ['Trunk', 'Branch', 'Private']:
            df = df.query(f'mutation == "{mutation_type}"')
            df = df.drop(columns='patient')
        elif mutation_type == 'Patient':
            # If mutation type is patient, we group by patient and
            # consider all mutated genes for each patient
            df = df.groupby('patient').agg({'genes': 'sum'}).reset_index()
            df = df.rename(columns={'patient': 'module'})
            df['genes'] = df['genes'].apply(lambda x: list(set(x)))

        df = df.explode('genes')

        modules[mutation_type] = df.groupby('genes')['module'].apply(set).to_dict()
    
    return (modules, cancer_type)

# Split the dataframe into chunks
print('Splitting the dataframe into chunks')
chunks = [group for cancer_type, group in df_complete.groupby('Cancertype')]

# Create a pool of workers
print('Creating a pool of workers')
pool = multiprocessing.Pool(num_processes)

# Compute the results
print('Computing the results')
results = pool.map(process_cancer_type, chunks)
    
# Close the pool
print('Closing the pool')
pool.close()

print('Constructing df')
# Make a dataframe with all modules for each gene and type of mutation e.g.
modules_df = pd.DataFrame()

for modules, cancer_type in results:
    for mutation_type, gene_modules in modules.items():
        for gene, modules in gene_modules.items():
            # Make tmp dataframe with gene, module, mutation type
            tmp = pd.DataFrame({'gene': gene, 'cancertype': cancer_type, 'module': list(modules), 'mutation': mutation_type})
            # Expand gene lists to one individual gene per row
            tmp = tmp.explode('gene')
            modules_df = pd.concat([modules_df, tmp])

# Correct gene names that were wrongly interpreted as dates
gene_symbols = {
	'1-Dec': 'DELEC1',
    '2-Mar': 'MARCHF2',
    '3-Mar': 'MARCHF3',
    '7-Mar': 'MARCHF7',
    'Mar-10': 'MARCHF10',
    '1-Mar': 'MARCHF1',
    '8-Mar': 'MARCHF8',
    '10-Mar': 'MARCHF10',
    '4-Mar': 'MARCHF4',
    '6-Mar': 'MARCHF6',
    'Mar-03': 'MARCHF3',
    'Mar-04': 'MARCHF4',
    '15-Sep': 'SEPTIN15',
    '4-Sep': 'SEPTIN4',
    '9-Sep': 'SEPTIN9',
    'Sep-08': 'SEPTIN8',
    '1-Sep': 'SEPTIN1',
    '3-Sep': 'SEPTIN3',
    '6-Sep': 'SEPTIN6',
    'Sep-09': 'SEPTIN9',
    'Sep-10': 'SEPTIN10',
    'Sep-11': 'SEPTIN11'
}

modules_df['gene'] = modules_df['gene'].replace(gene_symbols)

# Save modules_df
print(f'Saving modules_df to {output_dir}')
os.makedirs(output_dir, exist_ok=True)
modules_df.to_csv(os.path.join(output_dir, 'mutations_df.csv'), index=False)

print('Done: get_mutations.py')
