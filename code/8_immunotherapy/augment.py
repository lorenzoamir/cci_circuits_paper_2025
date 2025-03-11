import pandas as pd
import numpy as np
from imblearn.over_sampling import SMOTENC
import os
import sys

seed = 42
np.random.seed(seed)

train = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/train.csv', index_col=0)
test = pd.read_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/test.csv', index_col=0)

# drop patient columns, they don't make sense in the context of data augmentation
train = train.drop(columns=['dataset_id', 'patient_name', 'batch', 'treatment_when'])

categoricals = ['tissue', 'therapy_type']
genes = [col for col in train.columns if col not in categoricals + ['response_NR']]

print('Train:', train.shape)
print('Test:', test.shape)

# Data augmentation

###### SMOTENC
print('SMOTENC')
categorical_features = ['tissue', 'therapy_type']
X_aug, y_aug = SMOTENC(categorical_features=categorical_features, random_state=seed).fit_resample(train.drop(columns=['response_NR']), train['response_NR'])
print(X_aug.head())
print(X_aug.shape)

###### Random shuffling of genes

# create some new data by substituting a gene with one from another patient with the same tissue, therapy_type and response
# Copy data
X_new = X_aug.copy()
y_new = y_aug.copy()

# Group by tissue, therapy_type, and response
grouped = train.groupby(['tissue', 'therapy_type', 'response_NR'])

for (tissue, therapy_type, response), indices in grouped.groups.items():
    idx = X_new.index.intersection(indices)
    group = X_new.loc[idx, genes]

    if group.shape[0] < 2:
        continue

    # Select 10% of genes for shuffling
    to_shuffle = np.random.choice(group.columns, int(0.1 * group.shape[1]), replace=False)

    # Shuffle selected columns
    shuffled_values = group[to_shuffle].to_numpy()
    np.random.shuffle(shuffled_values)  # Shuffle in-place

    # Assign shuffled values back
    X_new.loc[idx, to_shuffle] = shuffled_values

# Concatenate efficiently
X_aug = pd.concat([X_aug, X_new], ignore_index=True)
y_aug = pd.concat([y_aug, y_new], ignore_index=True)
print(X_aug.shape)

###### Randomly drop features
X_new = X_aug.copy()
y_new = y_aug.copy()

# Set 10% of the cells in the new data to NaN (excluding the response column)
X_new = X_new.mask(np.random.random(X_new.shape) < 0.1)
# Add back the response column
X_new['response_NR'] = y_new
print(X_new.head())

X_aug = pd.concat([X_aug, X_new], ignore_index=True)
y_aug = pd.concat([y_aug, y_new], ignore_index=True)
print(X_aug.shape)

# Merge and save
X_aug['response_NR'] = y_aug
X_aug.to_csv('/home/lnemati/pathway_crosstalk/data/immunotherapy/train_aug.csv')

print('Done: augment.py')
