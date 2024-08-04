import pandas as pd
import numpy as np
import re
import os
import multiprocessing
from tqdm import tqdm

def extract_gene_symbols(mutation_string):
    if pd.isna(mutation_string) or mutation_string == '':
        return []

    # Split the string into individual entries
    gene_entries = mutation_string.split(',')

    # Filter entries to only those with parentheses
    filtered_entries = [entry for entry in gene_entries if '(' in entry and ')' in entry]

    # Remove the parentheses and any content within them
    cleaned_entries = [re.sub(r'\([^()]*\)', '', entry).strip() for entry in filtered_entries]

    # Return the cleaned entries
    return cleaned_entries

def jaccard_similarity(set1, set2):
    intersection = len(set1 & set2)
    union = len(set1 | set2)
    return intersection / union if union != 0 else 0

def calculate_jaccard_for_gene(gene1, all_genes, modules, progress_queue):
    similarities = np.zeros(len(all_genes))
    set1 = modules[gene1]
    for idx, gene2 in enumerate(all_genes):
        set2 = modules[gene2]
        similarities[idx] = jaccard_similarity(set1, set2)
    progress_queue.put(1)
    return similarities

def calculate_jaccard_matrix(all_genes, modules, num_processes):
    manager = multiprocessing.Manager()
    progress_queue = manager.Queue()
    total_tasks = len(all_genes)
    with multiprocessing.Pool(num_processes) as pool:
        results = pool.starmap_async(
            calculate_jaccard_for_gene,
            [(gene, all_genes, modules, progress_queue) for gene in all_genes]
        )

        with tqdm(total=total_tasks) as pbar:
            completed = 0
            while completed < total_tasks:
                progress_queue.get()
                completed += 1
                pbar.update(1)

        pool.close()
        pool.join()

    return np.array(results.get())

#gene_names = pd.read_csv('/home/lnemati/resources/biomart/ensembl_to_symbol_filtered.csv.gz')['Gene name'].values
output_dir = '/home/lnemati/pathway_crosstalk/results/coevolution/tumor'
df_complete = pd.read_csv('/projects/bioinformatics/DB/CancerTracer/intra_heterogeneity_data.txt', sep='\t')
df_complete = df_complete.set_index('Sample')
df_complete = df_complete.replace({'-': ''})

# Remove wrongly formatted row
df_complete = df_complete[~df_complete.index.isin(['p075_p3'])]

df_complete['Trunk'] = df_complete['Trunk_mutation'].apply(extract_gene_symbols)
df_complete['Branch'] = df_complete['Branch_mutation'].apply(extract_gene_symbols)
df_complete['Private'] = df_complete['Private_mutation'].apply(extract_gene_symbols)

# Only iterate on cancertypes with at least min_patients patients
min_patients = 10

mutation_types = ['All', 'Trunk', 'Branch', 'Private', 'Patient']

modules = {}

print('Running for each cancertype')
for cancer_type, group in df_complete.groupby('Cancertype'):

    print(f'Running for {cancer_type}')
    group = group[['Trunk', 'Branch', 'Private']]
    group = group.stack().reset_index()
    group = group.rename(columns={'level_1': 'mutation', 0: 'genes', 'Sample': 'patient'})
    group['module'] = group['patient'] + '_' + group['mutation']
    group = group[['genes', 'module', 'patient', 'mutation']]

    for mutation_type in mutation_types:
        df = group.copy()

        # Check if there are enough patients
        if len(df['patient'].unique()) < min_patients:
            print(f'Not enough patients for {cancer_type} {mutation_type}')
            continue

        if mutation_type in ['Trunk', 'Branch', 'Private']:
            df = df.query(f'mutation == "{mutation_type}"')
            df = df.drop(columns='patient')
        elif mutation_type == 'All':
            # If mutation type is all, we take all mutations,
            # but the co-occurrence has to be in the same category
            df = df.drop(columns='patient')
        elif mutation_type == 'Patient':
            # If mutation type is patient, we group by patient and
            # consider all mutated genes for each patient
            df = df.groupby('patient').agg({'genes': 'sum'}).reset_index()
            df = df.rename(columns={'patient': 'module'})
            df['genes'] = df['genes'].apply(lambda x: list(set(x)))

        df = df.explode('genes')

        modules[mutation_type] = df.groupby('genes')['module'].apply(set).to_dict()
        #modules[mutation_type] = {key: val for key, val in modules[mutation_type].items() if key in gene_names}

        all_genes = modules[mutation_type].keys()
        df = df[df['genes'].isin(all_genes)]

        num_processes = multiprocessing.cpu_count()

        print('Counting number of mutations for each gene')
        gene_counts = df['genes'].value_counts()

        print('Calculating jaccard matrix')
        jaccard_matrix = calculate_jaccard_matrix(all_genes, modules[mutation_type], num_processes)

        try:
            print('Making dataframe')
            jaccard_df = pd.DataFrame(jaccard_matrix, index=all_genes, columns=all_genes)

            # Saving the Jaccard index matrix and the number of mutations for each gene
            output_path = os.path.join(output_dir, cancer_type, mutation_type)
            os.makedirs(output_path, exist_ok=True)

            jaccard_df.to_csv(os.path.join(output_path, f'jaccard.csv'.lower()))
            gene_counts.to_csv(os.path.join(output_path, f'gene_counts.csv'.lower()))
        except Exception as e:
            print('Some error occurred in:', cancer_type, mutation_type)
            print('Perhaps there are no mutations for this combination')
            print('Error:')
            print(e)

# Run again but this time consider all cancertypes
print('Running for all cancertypes')
df_complete['Cancertype'] = 'All Cancers'

group = df_complete[['Trunk', 'Branch', 'Private']]
group = group.stack().reset_index()
group = group.rename(columns={'level_1': 'mutation', 0: 'genes', 'Sample': 'patient'})
group['module'] = group['patient'] + '_' + group['mutation']
group = group[['genes', 'module', 'patient', 'mutation']]

for mutation_type in mutation_types:
    df = group.copy()

    if mutation_type in ['Trunk', 'Branch', 'Private']:
        df = df.query(f'mutation == "{mutation_type}"')
        df = df.drop(columns='patient')
    elif mutation_type == 'All':
        # If mutation type is all, we take all mutations,
        # but the co-occurrence has to be in the same category
        df = df.drop(columns='patient')
    elif mutation_type == 'Patient':
        # If mutation type is patient, we group by patient and
        # consider all mutated genes for each patient
        df = df.groupby('patient').agg({'genes': 'sum'}).reset_index()
        df = df.rename(columns={'patient': 'module'})
        df['genes'] = df['genes'].apply(lambda x: list(set(x)))

    df = df.explode('genes')

    modules[mutation_type] = df.groupby('genes')['module'].apply(set).to_dict()
    #modules[mutation_type] = {key: val for key, val in modules[mutation_type].items() if key in gene_names}

    all_genes = modules[mutation_type].keys()
    df = df[df['genes'].isin(all_genes)]

    print('Counting number of mutations for each gene')
    gene_counts = df['genes'].value_counts()

    num_processes = multiprocessing.cpu_count()

    print('Calculating matrices')
    # Make a matrix with the Jaccard index
    # Also make a matrix with the number of co-occurrences (intersection of the sets of mutations)
    jaccard_matrix = calculate_jaccard_matrix(all_genes, modules[mutation_type], num_processes)

    print('Making dataframe')
    jaccard_df = pd.DataFrame(jaccard_matrix, index=all_genes, columns=all_genes)

    

    # Saving the Jaccard index matrix and the number of mutations for each gene
    output_path = os.path.join(output_dir, 'All Cancers', mutation_type)
    os.makedirs(output_path, exist_ok=True)

    jaccard_df.to_csv(os.path.join(output_path, f'jaccard.csv'.lower()))
    gene_counts.to_csv(os.path.join(output_path, f'gene_counts.csv'.lower()))

print("Done: coevolution_jaccard.py")
