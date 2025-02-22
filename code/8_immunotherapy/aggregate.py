import os
import pandas as pd

# Read all .csv files in the data_dir and merge them into a single .csv file
data_dir = '/home/lnemati/pathway_crosstalk/results/immunotherapy/aurocs'
result_dir = '/home/lnemati/pathway_crosstalk/results/immunotherapy'

dfs = []
for file in os.listdir(data_dir):
    if file.endswith('.csv'):
        df = pd.read_csv(os.path.join(data_dir, file))
        dfs.append(df)

df = pd.concat(dfs, ignore_index=True)

df.to_csv(os.path.join(result_dir, 'all_aurocs.csv'), index=False)

print('Done: aggregate.py')
