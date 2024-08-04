import pandas as pd
import multiprocessing as mp
import numpy as np
import os

output_path = '/home/lnemati/pathway_crosstalk/data/tumor_coev'

# Read mutations df
print('Reading mutations df')
df = pd.read_csv('/home/lnemati/pathway_crosstalk/results/coevolution/tumor/mutations_df.csv')

print(df.head())

min_modules = 30

# Define a function to apply to each group
def save_group(group):
    (cancertype, mutationtype), data = group
    
    # For each group, determine if the cancertype has enough modules
    n_modules = df[df['cancertype'] == cancertype]['module'].nunique()
    
    if n_modules < min_modules:
        return

    # If not exist create tissue and mutationtype folders
    os.makedirs(os.path.join(output_path, cancertype), exist_ok=True)
    os.makedirs(os.path.join(output_path, cancertype, mutationtype), exist_ok=True) 

    data.to_csv(f'{output_path}/{cancertype}/{mutationtype}/mutations_df.csv', index=False)

# Function to perform multiprocessing with groupby
def parallel_groupby(df):
    # Split the DataFrame into groups by multiple columns
    groups = df.groupby(['cancertype', 'mutation'])

    # Get number of cpus from NCPUS environment variable
    ncpus = int(os.getenv('NCPUS', mp.cpu_count()))

    # Create a pool of worker processes
    with mp.Pool(processes=ncpus) as pool:
        # Save all groups
        pool.map(save_group, groups)


df = df.drop_duplicates()

parallel_groupby(df)

print('Done: separate.py')
