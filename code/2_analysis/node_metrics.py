import pandas as pd
import numpy as np
import gseapy as gp
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Compare hub genes from WGCNAs')

parser.add_argument('--dir-list', type=str, help='List of directories to compare, delimited by commas')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results/hubs'
# Create output directory if it does not exist
os.makedirs(output, exist_ok=True)

perc = 0.05

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]

check_files = [
    'degree_df.csv.gz',
    'interactions.csv',
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

all_genes = set()
all_interactors = set()

normal_tissues = []
tumor_tissues = []

metrics = ['degree', 'closeness', 'betweenness', 'pagerank', 'clustering_local']
top_genes = {metric: {} for metric in metrics}

print('Reading all files')
for d in dir_list:

    print(d)
    name, condition = get_name_condition(d)
    print(name)
   
    if condition == 'normal':
        normal_tissues.append(name)
    elif condition == 'tumor':
        tumor_tissues.append(name)
   
    degree_df = pd.read_csv(os.path.join(d, 'degree_df.csv.gz'), index_col=0)
    interactions = pd.read_csv(os.path.join(d, 'interactions.csv'))

    # Add all genes to set
    all_genes.update(set(degree_df.index))
    
    # Create set of all interactors by taking all unique values that appear in 
    # columns interactor1, interactor2, ..., interactor7
    for i in range(1, 8):
        all_interactors.update(set(interactions['interactor{}'.format(i)]))

    print('Number of interactors: {}'.format(len(all_interactors)))

    for metric in metrics:
        # Get top genes
        print('Getting top genes by {}'.format(metric))
        #hubs[name] = degree_df.sort_values('degree', ascending=False).head(int(perc * len(degree_df))).index
        top_genes[metric][name] = degree_df.sort_values(metric, ascending=False).head(int(perc * len(degree_df))).index

print('Done reading all files')

print('Scoring each gene')
# Score each gene by how many times it appears in the top 5% of genes
columns = ['f_{}'.format(metric) + '_tumor' for metric in metrics] + \
          ['f_{}'.format(metric) + '_normal' for metric in metrics] + \
          ['interactor']
genes_scores = pd.DataFrame(index=list(all_genes), columns=columns)

for gene in all_genes:
    # Mark CCC genes
    genes_scores.loc[gene, 'interactor'] = 1 if gene in all_interactors else 0

    for metric in metrics:
        genes_scores.loc[gene, 'f_{}'.format(metric) + '_tumor'] = sum(
            [gene in top_genes[metric][name] for name in tumor_tissues]
        ) / len(tumor_tissues)
        genes_scores.loc[gene, 'f_{}'.format(metric) + '_normal'] = sum(
            [gene in top_genes[metric][name] for name in normal_tissues]
        ) / len(normal_tissues) 

print('Done scoring each gene')
   
# Save results
genes_scores.to_csv(os.path.join(output, 'n_tissues_top_node_metrics.csv'))

# ----- Enrichment analysis -----

reactome = '/home/lnemati/resources/reactome/ReactomePathways.gmt'
hallmarks = '/home/lnemati/resources/cancer_hallmarks/gsea_hallmarks.gmt'

# Make previous code into a function
def enrichment_analysis(genes_scores, metric, output, genesets, subset=False):
    # Create directory if it does not exist and save results
    os.makedirs(os.path.join(output, metric), exist_ok=True)

    # Get genes with higher score in tumor than in normal
    tumor_genes = genes_scores[
        (genes_scores['f_{}'.format(metric) + '_tumor'] > genes_scores['f_{}'.format(metric) + '_normal'])
    ].index.tolist()

    try:
        enr = gp.enrichr(
            gene_list=tumor_genes,
            gene_sets=genesets,
            background=list(all_genes),
            outdir=None,
        )

        enr = enr.results
        if subset:
            enr = enr[enr["Adjusted P-value"] < 0.05]
        enr.to_csv(os.path.join(output, metric, 'enrichment_tumor.csv'), index=False)

    except:
        print(f"Error in enrichment analysis for tumor {metric}", file=sys.stderr)
        print(f"Error in enrichment analysis for tumor {metric}", file=sys.stdout)

    # Repeat for normal
    normal_genes = genes_scores[
        (genes_scores['f_{}'.format(metric) + '_tumor'] < genes_scores['f_{}'.format(metric) + '_normal'])
    ].index.tolist()

    try:
        enr = gp.enrichr(
            gene_list=normal_genes,
            gene_sets=genesets,
            background=list(all_genes),
            outdir=None,
        )
    except:
        print(f"Error in enrichment analysis for normal {metric}", file=sys.stderr)
        print(f"Error in enrichment analysis for normal {metric}", file=sys.stdout)

    enr = enr.results
    if subset:
        enr = enr[enr["Adjusted P-value"] < 0.05]

    # Create directory if it does not exist and save results
    os.makedirs(os.path.join(output, metric), exist_ok=True)
    enr.to_csv(os.path.join(output, metric, 'enrichment_normal.csv'), index=False)

    print(f'Finished enrichment analysis for {metric}')

# Run reactome enrichment analysis for each metric
all_hubs_output = os.path.join(output, 'enrichment', 'reactome', 'all_hubs')
for metric in metrics:
    enrichment_analysis(genes_scores, metric, all_hubs_output, reactome, subset=True)
# Subset to CCC (interactor) genes and rerun enrichment analysis
ccc_output = os.path.join(output, 'enrichment', 'reactome', 'ccc_only')
ccc_genes_scores = genes_scores[genes_scores['interactor'] == 1]
for metric in metrics:
    enrichment_analysis(ccc_genes_scores, metric, ccc_output, reactome, subset=True)

# Run hallmarks enrichment analysis for each metric
all_hubs_output = os.path.join(output, 'enrichment', 'hallmarks', 'all_hubs')
for metric in metrics:
    enrichment_analysis(genes_scores, metric, all_hubs_output, hallmarks, subset=False)
# Subset to CCC (interactor) genes and rerun enrichment analysis
ccc_output = os.path.join(output, 'enrichment', 'hallmarks', 'ccc_only')
ccc_genes_scores = genes_scores[genes_scores['interactor'] == 1]
for metric in metrics:
    enrichment_analysis(ccc_genes_scores, metric, ccc_output, hallmarks, subset=False)

print('Done: node_metrics.py')
