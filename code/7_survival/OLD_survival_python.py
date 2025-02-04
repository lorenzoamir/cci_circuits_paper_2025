# Import necessary libraries
import pandas as pd
import numpy as np
import lifelines
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
import os
import sys
import re
import argparse

# Define constants
MIN_PATIENTS = 10

# Get number of CPUs
ncpus = os.environ.get('NCPUS')

# Parse arguments
parser = argparse.ArgumentParser(description='Survival analysis for all patients')

parser.add_argument('--tissuefile', type=str, help='Path to tissue file')
parser.add_argument('--output_dir', type=str, help='Path to output directory')

args = parser.parse_args()

# Read arguments
tissuefile = args.tissuefile
output_dir = args.output_dir

# Read the survival data
print('Reading survival data')
df = pd.read_csv(tissuefile, index_col=0)

# Numbers after the last - in the patient ids represent samples, average all samples for each patient
df['patient'] = [x.rsplit('-', 1)[0] for x in df.index]

# condition,tissue,type,gender,study,OS,OS.time,DSS,DSS.time,DFI,DFI.time,PFI,PFI.time
# some columns are not expression values, don't average them
aggregates = {k: 'first' for k in ['condition', 'tissue', 'type', 'gender', 'study', 'OS', 'OS.time', 'DSS', 'DSS.time', 'DFI', 'DFI.time', 'PFI', 'PFI.time']}
patient_info = df.groupby('patient').agg(aggregates)

# Other columns are expression values, average them
expression = df.drop(columns=aggregates.keys())
expression = expression.groupby('patient').mean()

# Combine the two dataframes
df = pd.concat([patient_info, expression], axis=1)

# Print number of patients (rows) in the dataset
print(f'Number of patients: {df.shape[0]}')

# Read motifs file and filter for specific types
motifs_file = '/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv'
motifs = pd.read_csv(motifs_file, index_col='Interaction')

motifs = motifs.index
print(f'Number of cliques (subset): {len(motifs)}')

# Get ccc interactions from motifs
ccc = motifs.str.split('&', expand=True)
ccc = set(ccc.get_level_values(0)).union(set(ccc.get_level_values(1)))
ccc = list(ccc)
print('Number of ccc interactions:', len(ccc))

def survival_analysis(interaction, df):
    print(interaction, file=sys.stdout)
    print(interaction, file=sys.stderr)

    # Split the interaction string into individual genes
    genes = list(set(re.split(r'[+&_]', interaction)))

    # Check if multiple conditions are present
    multiple = 'condition' in df.columns and len(df['condition'].unique()) > 1

    # Select relevant columns
    cols = ['OS.time', 'OS', 'condition'] + genes if multiple else ['OS.time', 'OS'] + genes
    df = df[cols].copy()

    # Convert OS time to years
    df['OS.time'] = df['OS.time'] / 365.0

    # Prepare dataframes for high and low expression groups
    high_expression_group = pd.DataFrame()
    low_expression_group = pd.DataFrame()

    print('Splitting')
    if multiple:
        # If multiple conditions are present, split by condition
        for condition in df['condition'].unique():
            tissue_df = df[df['condition'] == condition].copy()
            medians = tissue_df[genes].median(axis=0, skipna=True)
            
            # Identify patients above and below the median for all genes
            above_all_genes = tissue_df[genes].gt(medians, axis=1).all(axis=1)
            below_all_genes = tissue_df[genes].le(medians, axis=1).all(axis=1)
            
            high_expression_group = pd.concat([high_expression_group, tissue_df[above_all_genes]])
            low_expression_group = pd.concat([low_expression_group, tissue_df[below_all_genes]])
    else:
        # If no condition column, split directly
        medians = df[genes].median(axis=0, skipna=True)
        above_all_genes = df[genes].gt(medians, axis=1).all(axis=1)
        below_all_genes = df[genes].le(medians, axis=1).all(axis=1)
        
        high_expression_group = df[above_all_genes]
        low_expression_group = df[below_all_genes]

    # Add a group column to indicate high (1) or low (0) expression
    high_expression_group['group'] = 1
    low_expression_group['group'] = 0

    # Combine groups
    df = pd.concat([high_expression_group, low_expression_group])

    print('Subsetting')
    df = df[df['group'].isin([0, 1])]

    # Count the number of patients in each group
    n_above_all = len(high_expression_group)
    n_below_all = len(low_expression_group)

    # If there are too few patients in either group, return default NA values
    if n_above_all < MIN_PATIENTS or n_below_all < MIN_PATIENTS:
        print('Too few patients')
        return {
            'hr': np.nan,
            'n_patients_low': n_below_all,
            'n_patients_high': n_above_all,
            'loglik_pval': np.nan,
            'logrank_pval': np.nan,
            'concordance_index': np.nan,
            'ci_low': np.nan,
            'ci_high': np.nan,
            'se': np.nan
        }

    print('Converting')
    # Select relevant columns for Cox model
    if multiple:
        df = df[['OS.time', 'OS', 'group', 'condition']]
    else:
        df = df[['OS.time', 'OS', 'group']]

    # Drop any rows with missing values
    df = df.dropna()

    # Fit Cox Proportional Hazards Model
    try:
        print('Fitting model')
        cph = CoxPHFitter(penalizer=0.1)
        if multiple:
            cph.fit(df, duration_col='OS.time', event_col='OS', strata='condition')
        else:
            cph.fit(df, duration_col='OS.time', event_col='OS', step_size=0.1)

        # Extract model metrics
        hr = np.exp(cph.hazard_ratios_['group'])
        loglik_pval = cph.log_likelihood_ratio_test().p_value
        logrank_pval = lifelines.statistics.logrank_test(
            df[df['group'] == 0]['OS.time'],
            df[df['group'] == 1]['OS.time'],
            event_observed_A=df[df['group'] == 0]['OS'],
            event_observed_B=df[df['group'] == 1]['OS']
        ).p_value
        concordance = concordance_index(df['OS.time'], -df['group'].astype(float), df['OS'])
        ci = cph.confidence_intervals_.loc['group']
        ci_low = np.exp(ci['95% lower-bound'])
        ci_high = np.exp(ci['95% upper-bound'])
        se = hr * cph.standard_errors_['group']

        return {
            'hr': hr,
            'n_patients_low': n_below_all,
            'n_patients_high': n_above_all,
            'logrank_pval': logrank_pval,
            'loglik_pval': loglik_pval,
            'concordance_index': concordance,
            'ci_low': ci_low,
            'ci_high': ci_high,
            'se': se
        }

    except Exception as e:
        print(f"Error during Cox model fitting: {e}")
        return {
            'hr': np.nan,
            'n_patients_low': n_below_all,
            'n_patients_high': n_above_all,
            'logrank_pval': np.nan,
            'loglik_pval': np.nan,
            'concordance_index': np.nan,
            'ci_low': np.nan,
            'ci_high': np.nan,
            'se': np.nan
        }

# Test all ccc interactions
ccc_results = []
for interaction in ccc:
    result = survival_analysis(interaction, df)
    # Append results with interaction details
    ccc_results.append({
        'interaction': interaction,
        'hr': result['hr'],
        'n_patients_low': result['n_patients_low'],
        'n_patients_high': result['n_patients_high'],
        'logrank_pval': result['logrank_pval'],
        'loglik_pval': result['loglik_pval'],
        'concordance_index': result['concordance_index'],
        'ci_low': result['ci_low'],
        'ci_high': result['ci_high'],
        'se': result['se'],
        'type': 'ccc'
    })

# Convert results to a DataFrame
ccc_results_df = pd.DataFrame(ccc_results)

# Test all crosstalk interactions
crosstalk_results = []
for interaction in motifs:
    result = survival_analysis(interaction, df)
    # Append results with interaction details
    crosstalk_results.append({
        'interaction': interaction,
        'hr': result['hr'],
        'n_patients_low': result['n_patients_low'],
        'n_patients_high': result['n_patients_high'],
        'logrank_pval': result['logrank_pval'],
        'loglik_pval': result['loglik_pval'],
        'concordance_index': result['concordance_index'],
        'ci_low': result['ci_low'],
        'ci_high': result['ci_high'],
        'se': result['se'],
        'type': 'crosstalk'
    })

# Convert results to a DataFrame
crosstalk_results_df = pd.DataFrame(crosstalk_results)

# Combine results
all_results = pd.concat([ccc_results_df, crosstalk_results_df], ignore_index=True)
all_results.sort_values(by='logrank_pval', inplace=True)

# Save results to a CSV file
tissue_name = os.path.basename(tissuefile).replace('.csv', '')
os.makedirs(output_dir, exist_ok=True)
output_file = os.path.join(output_dir, f'{tissue_name}.csv')

print(f'Saving results to: {output_file}')
all_results.to_csv(output_file, index=False)

print('Done: survival_analysis.py')
