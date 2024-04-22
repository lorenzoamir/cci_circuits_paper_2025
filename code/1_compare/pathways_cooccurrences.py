import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Check all dataframes to find pathway co-occurrences')

parser.add_argument('-o', '--outputdir', type=str, help='Path to output file')

args = parser.parse_args()

if args.outputdir is None:
    args.outputdir = '/home/lnemati/pathway_crosstalk/results/'

# Searching for pathways_cooccurrences.csv.gz files
parent_dir = '/projects/bioinformatics/DB/Xena/TCGA_GTEX/by_tissue_primary_vs_normal'

t_pw_files = []
n_pw_files = []

t_enr_files = []
n_enr_files = []

print('Searching for pathways_cooccurrences.csv files')
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.startswith('pathways_cooccurrences.csv'):
            if '/tumor' in root:
                t_pw_files.append(os.path.join(root, file))
            elif '/normal' in root:
                n_pw_files.append(os.path.join(root, file))

print('Searching for module_enrichment.csv files')
for root, dirs, files in os.walk(parent_dir):
    for file in files:
        if file.startswith('module_enrichment.csv'):
            if '/tumor' in root:
                t_enr_files.append(os.path.join(root, file))
            elif '/normal' in root:
                n_enr_files.append(os.path.join(root, file))

def get_all_indices(paths):
    # Read one file at a time and get all indices
    all_indices = set()
    for path in paths:
        df = pd.read_csv(path, index_col='Term')
        all_indices.update(df.index)

    return list(all_indices)

all_indices = get_all_indices(t_enr_files + n_enr_files)
print('All pathways: ', len(all_indices))

print('First 5 pathways: ')
print(all_indices[:5])
print()

# ----- Genes overlap -----
print('Genes overlap')

# Calculate pathways gene sets jaccaard similarity only once and then
# use the same values for both tumor and normal

pw_genes = {}
for key in all_indices:
    pw_genes[key] = set()

for path in t_enr_files + n_enr_files:
    print('Reading file: ', path)
    df = pd.read_csv(path, index_col='Term')[['Genes']]
    for index, row in df.iterrows():
        print(index)
        pw_genes[index].update(row['Genes'].split(';'))

print('First 3 pathways genes: ')

for key in list(pw_genes.keys())[:3]:
    print(key)
    print(pw_genes[key])

genes_jaccard_df = pd.DataFrame(0, index=all_indices, columns=all_indices)

def jaccard_similarity(set1, set2):
    if len(set1) == 0 or len(set2) == 0:
        return 0
    return len(set1.intersection(set2)) / len(set1.union(set2))

for i in range(len(all_indices)):
    for j in range(i+1, len(all_indices)):
        index1 = all_indices[i]
        index2 = all_indices[j]

        genes_jaccard_df.at[index1, index2] = jaccard_similarity(pw_genes[index1], pw_genes[index2])
        genes_jaccard_df.at[index2, index1] = genes_jaccard_df.at[index1, index2]

genes_jaccard_df.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'pathways_genesets_jaccard.csv'))
print()

# ----- Tumor -----
print('Tumor')

print('Counting co-occurrences for tumor files')
t_sum_df = pd.DataFrame(0, index=all_indices, columns=all_indices)

for path in t_pw_files:
    df = pd.read_csv(path, index_col=0)
    t_sum_df = t_sum_df.add(df, fill_value=0)

# Create pathways_cooccurrences directory if it doesn't exist
if not os.path.exists(os.path.join(args.outputdir, 'pathways_cooccurrences')):
    os.makedirs(os.path.join(args.outputdir, 'pathways_cooccurrences'))

t_sum_df.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'tumor_cooccurrences_counts.csv'))

print('Calculating Jaccard similarity')

# For each pathway in all_indices, get the set of modules that are enriched in that pathway
# Then calculate the Jaccard similarity between the set of modules for each pair of pathways

t_pw_modules = {}

for path in t_enr_files:
    name = path.split('/')[-2]
    df = pd.read_csv(path, index_col='Term')[['module']]
    for index, row in df.iterrows():
        if index not in t_pw_modules:
            t_pw_modules[index] = set()
        # Add name to the module to differentiate between different tissues
        t_pw_modules[index].add(name + '_' + str(row['module']))

print('First 3 pathways modules: ')
for key in list(t_pw_modules.keys())[:3]:
    print(key)
    print(t_pw_modules[key])

t_modules_jaccard_df = pd.DataFrame(0, index=all_indices, columns=all_indices)

for i in range(len(all_indices)):
    for j in range(i+1, len(all_indices)):
        index1 = all_indices[i]
        index2 = all_indices[j]

        # if one of the pathways is not in the dictionary, set the Jaccard similarity to 0
        if index1 not in t_pw_modules.keys() or index2 not in t_pw_modules.keys():
            t_modules_jaccard_df.at[index1, index2] = 0.
            t_modules_jaccard_df.at[index2, index1] = 0.
        else:
            t_modules_jaccard_df.at[index1, index2] = jaccard_similarity(t_pw_modules[index1], t_pw_modules[index2])
            t_modules_jaccard_df.at[index2, index1] = t_modules_jaccard_df.at[index1, index2]

t_modules_jaccard_df.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'tumor_modules_jaccard.csv'))
print()

# ----- Normal -----
print('Normal')

print('Counting co-occurrences for normal files')
n_sum_df = pd.DataFrame(0, index=all_indices, columns=all_indices)

for path in n_pw_files:
    df = pd.read_csv(path, index_col=0)
    n_sum_df = n_sum_df.add(df, fill_value=0)

n_sum_df.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'normal_cooccurrences_counts.csv'))

print('Calculating Jaccard similarity')

n_pw_modules = {}

for path in n_enr_files:
    name = path.split('/')[-3]
    df = pd.read_csv(path, index_col='Term')[['module']]
    for index, row in df.iterrows():
        if index not in n_pw_modules:
            n_pw_modules[index] = set()
        n_pw_modules[index].add(name + '_' + str(row['module']))

print('First 3 pathways modules: ')

for key in list(n_pw_modules.keys())[:3]:
    print(key)
    print(n_pw_modules[key])

n_modules_jaccard_df = pd.DataFrame(0, index=all_indices, columns=all_indices)

for i in range(len(all_indices)):
    for j in range(i+1, len(all_indices)):
        index1 = all_indices[i]
        index2 = all_indices[j]

        if index1 not in n_pw_modules.keys() or index2 not in n_pw_modules.keys():
            n_modules_jaccard_df.at[index1, index2] = 0.
            n_modules_jaccard_df.at[index2, index1] = 0.
        else:
            n_modules_jaccard_df.at[index1, index2] = jaccard_similarity(n_pw_modules[index1], n_pw_modules[index2])
            n_modules_jaccard_df.at[index2, index1] = n_modules_jaccard_df.at[index1, index2]

n_modules_jaccard_df.to_csv(os.path.join(args.outputdir, 'pathways_cooccurrences', 'normal_modules_jaccard.csv'))
print()

print('Done: pathways_cooccurrences.py')
