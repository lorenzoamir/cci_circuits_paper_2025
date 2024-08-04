import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import os
import sys
import argparse
import multiprocessing

parser = argparse.ArgumentParser(description='Compare coexpression and coevolution')

parser.add_argument('--dir-list', type=str, help='List of directories', required=True)
parser.add_argument('--use-existing', action='store_true', help='Use existing coexpression and coevolution matrices')

args = parser.parse_args()

dir_list = args.dir_list

output = '/home/lnemati/pathway_crosstalk/results/coevolution'
if not os.path.exists(output):
    os.makedirs(output)

logfile = './coevolution.log'

dir_list = dir_list.split(',')
dir_list = [x.strip() for x in dir_list]

check_files = [
    'tom.csv.gz',
]

def printall(string):
    print(string, file=sys.stdout)
    print(string, file=sys.stderr)
    print(string, file=open(logfile, 'a'))

# Sanity check
for d in dir_list:
    if not os.path.exists(d):
        # error message to standard error and standard output
        printall('Error: Directory {} does not exist'.format(d))
    # Check if all files are present
    for f in check_files:
        if not os.path.exists(os.path.join(d, f)):
            printall('Error: File {} does not exist in directory {}'.format(f, d))

printall('Reading coevolution matrix')
coev = pd.read_csv('/home/lnemati/resources/coevolution/jaccard_genes.csv.gz', index_col=0)
printall('Converting coevolution matrix to float32')
coev = coev.astype(np.float32)
printall(coev.shape)    

def get_name_condition(d):
    if '/tumor/' in d: 
        condition = 'tumor'
    elif '/normal/' in d:
        condition = 'normal'
    else:
        printall('Error: Directory {} does not contain /tumor or /normal'.format(d))
        return None, None
    
    name = d.split('/')[-1]
    return name, condition

# Split dir_list into normal and tumor dirs
normal_dirs = [d for d in dir_list if '/normal/' in d]
tumor_dirs = [d for d in dir_list if '/tumor/' in d]

if not args.use_existing:
    for cond in ['normal', 'tumor']:
        printall('Condition: {}'.format(cond))
        if cond == 'normal':
            printall('Normal')
            cond_dir_list = normal_dirs
            output_tensor = np.zeros((coev.shape[0], coev.shape[1], len(normal_dirs)), dtype=np.float32)
            output_tensor[:] = np.nan
        elif cond == 'tumor':
            printall('Tumor')
            cond_dir_list = tumor_dirs
            output_tensor = np.zeros((coev.shape[0], coev.shape[1], len(tumor_dirs)), dtype=np.float32)
            output_tensor[:] = np.nan

        # names and conditions
        names = []

        for i, d in enumerate(cond_dir_list):
            printall('Directory {}'.format(d))
            name, condition = get_name_condition(d)
        
            printall('Name: {}'.format(name))
            printall('Condition: {}'.format(condition))
            if name is None:
                continue
            names.append(name)

            # Read the TOM, index are gene names (strings)
            printall('Reading TOM')
            df = pd.read_csv(os.path.join(d, 'tom.csv.gz'), index_col=0)
            df = df.astype(np.float32)
            
            # Expand tom to the full size of the coevolution matrix
            printall('Expanding TOM')
            common_genes = coev.index.intersection(df.index)
            all_genes = coev.index
            tom = pd.DataFrame(np.nan, index=all_genes, columns=all_genes)
            tom.loc[common_genes, common_genes] = df.loc[common_genes, common_genes]

            # Fill current slice of the output tensor with values from the TOM
            printall(f'Filling tensor slice {i}')
            output_tensor[:, :, i] = tom.values

        # output path
        printall('Saving')
        tensor_path = os.path.join(output, 'coexpression_tensor_{}.npz'.format(cond))

        # Save the tensor
        printall('Saving tensor')
        np.savez_compressed(tensor_path, tensor=output_tensor, names=names)

        # Collapse the tensor by averaging over each condition, discard nans
        printall('Collapsing tensor')
        mean = np.nanmean(output_tensor, axis=2)
        std = np.nanstd(output_tensor, axis=2)
        maxval = np.nanmax(output_tensor, axis=2)
        
        # Save the collapsed tensors as a dataframes
        printall('Saving collapsed tensor')
        mean = pd.DataFrame(mean, index=coev.index, columns=coev.columns, dtype=np.float32)
        std = pd.DataFrame(std, index=coev.index, columns=coev.columns, dtype=np.float32)
        maxval = pd.DataFrame(maxval, index=coev.index, columns=coev.columns, dtype=np.float32)

        mean.to_csv(os.path.join(output, 'avg_coexpression_{}.csv.gz'.format(cond)))
        std.to_csv(os.path.join(output, 'std_coexpression_{}.csv.gz'.format(cond)))
        maxval.to_csv(os.path.join(output, 'max_coexpression_{}.csv.gz'.format(cond)))

# Function to calculate correlations for a single gene
def correlation(gene):
    try:
        coev_gene = coev.loc[gene]

        # Get the coexpression values for the gene
        mean_gene = mean.loc[gene]
        std_gene = std.loc[gene]
        max_gene = maxval.loc[gene]

        # Calculate correlations
        corr_mean, pval_mean = spearmanr(coev_gene, mean_gene, nan_policy='omit')
        corr_std, pval_std = spearmanr(coev_gene, std_gene, nan_policy='omit')
        corr_max, pval_max = spearmanr(coev_gene, max_gene, nan_policy='omit')

        return gene, corr_mean, pval_mean, corr_std, pval_std, corr_max, pval_max
    except Exception as e:
        print(f"Error processing gene {gene}: {e}")
        return gene, None, None, None, None, None, None
for cond in ['normal', 'tumor']:
    printall('Condition: {}'.format(cond))
    # Create dataframe with genes as index
    df = pd.DataFrame(index=coev.index, columns=['corr_mean', 'corr_mean_pval', 'corr_std', 'corr_std_pval', 'corr_max', 'corr_max_pval'])

for cond in ['normal', 'tumor']:
    printall('Condition: {}'.format(cond))
    # Load the mean, std and max coexpression matrices
    mean = pd.read_csv(os.path.join(output, 'avg_coexpression_{}.csv.gz'.format(cond)), index_col=0)
    std = pd.read_csv(os.path.join(output, 'std_coexpression_{}.csv.gz'.format(cond)), index_col=0)
    maxval = pd.read_csv(os.path.join(output, 'max_coexpression_{}.csv.gz'.format(cond)), index_col=0)

    # Create a pool of workers
    with multiprocessing.Pool() as pool:
        results = pool.map(correlation, coev.index)

    # Assign the results back to the DataFrame
    for gene, corr_mean, pval_mean, corr_std, pval_std, corr_max, pval_max in results:
        if corr_mean is not None:  # Check if the result is valid
            df.loc[gene, 'corr_mean'] = corr_mean
            df.loc[gene, 'corr_mean_pval'] = pval_mean
            df.loc[gene, 'corr_std'] = corr_std
            df.loc[gene, 'corr_std_pval'] = pval_std
            df.loc[gene, 'corr_max'] = corr_max
            df.loc[gene, 'corr_max_pval'] = pval_max
        else:
            print(f"Skipping gene {gene} due to an error.")

    # Save the dataframe
    df.to_csv(os.path.join(output, 'coexpression_correlation_{}.csv.gz'.format(cond)))

printall('Done: coevolution.py')

## Using a context manager for the multiprocessing pool
#with multiprocessing.Pool() as pool:
#    results = pool.map(correlation, coev.index)
#
## Assign the results back to the DataFrame
#for gene, corr_mean, pval_mean, corr_std, pval_std, corr_max, pval_max in results:
#    if corr_mean is not None:  # Check if the result is valid
#        df.loc[gene, 'corr_mean'] = corr_mean
#        df.loc[gene, 'corr_mean_pval'] = pval_mean
#        df.loc[gene, 'corr_std'] = corr_std
#        df.loc[gene, 'corr_std_pval'] = pval_std
#        df.loc[gene, 'corr_max'] = corr_max
#        df.loc[gene, 'corr_max_pval'] = pval_max
#    else:
#        print(f"Skipping gene {gene} due to an error.")

    ## For each gene make correlation with coevolution
    #for gene in coev.index:
    #    coev_gene = coev.loc[gene]

    #    # Get the coexpression values for the gene
    #    mean_gene = mean.loc[gene]
    #    std_gene = std.loc[gene]
    #    max_gene = maxval.loc[gene]

    #    # Calculate correlation
    #    corr_mean, pval_mean = spearmanr(coev_gene, mean_gene, nan_policy='omit')
    #    corr_std, pval_std = spearmanr(coev_gene, std_gene, nan_policy='omit')
    #    corr_max, pval_max = spearmanr(coev_gene, max_gene, nan_policy='omit')

    #    # Save the values
    #    df.loc[gene, 'corr_mean'] = corr_mean
    #    df.loc[gene, 'corr_mean_pval'] = pval_mean
    #    df.loc[gene, 'corr_std'] = corr_std
    #    df.loc[gene, 'corr_std_pval'] = pval_std
    #    df.loc[gene, 'corr_max'] = corr_max
    #    df.loc[gene, 'corr_max_pval'] = pval_max

# Save the dataframe
#df.to_csv(os.path.join(output, 'coexpression_correlation_{}.csv.gz'.format(cond)))
