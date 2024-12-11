import pandas as pd
import os
import sys
import numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index
from sklearn.decomposition import PCA
import re
import argparse

# Parse arguments
parser = argparse.ArgumentParser(description='Survival analysis for all patients')

parser.add_argument('--tissuefile', type=str, help='Path to tissue file')

args = parser.parse_args()

# Read survival data
print('Reading survival data')
df = pd.read_csv(args.tissuefile, index_col=0)
print('Number of patients:', df.shape[0])

# Read motifs
motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv', index_col='Interaction')
motifs = motifs.query('Type in ["3_clique", "4_clique", "4_no_crosstalk"]').index
print('Number of cliques:', len(motifs))

# Get ccc interactions from motifs
ccc = motifs.str.split('&', expand=True)
ccc = set(ccc.get_level_values(0)).union(set(ccc.get_level_values(1)))
ccc = list(ccc)
print('Number of ccc interactions:', len(ccc))

def survival(interactions, df):  
    print(interactions, file=sys.stdout)
    print(interactions, file=sys.stderr)
    # Extract genes from the interaction string
    genes = list(set(re.split(r'[+&_]', interactions)))
    
    # Check if multiple tissues are present
    pan_cancer = True if len(df['tissue'].unique()) > 1 else False

    # Filter and process the DataFrame
    cols = ['OS.time', 'OS', 'tissue'] + genes if pan_cancer else ['OS.time', 'OS'] + genes

    df = df[cols]
    # Convert columns to numeric, except for tissue
    df.loc[:, [col for col in df.columns if col != 'tissue']] = df.loc[:, [col for col in df.columns if col != 'tissue']].apply(pd.to_numeric, errors='coerce')
    df = df.dropna()

    # Convert OS time to years
    df['OS.time'] = df['OS.time'] / 365
    
    if not pan_cancer:
        # If there is only one tissue, split the patients into high and low expression groups
        median = df[genes].quantile(0.5)
        high_expression_group = df[(df[genes] > median).T.all()]
        low_expression_group = df[(df[genes] <= median).T.all()] 
    else:
        # If there are multiple tissues, split the patients into high and low expression groups for each tissue
        high_expression_group = pd.DataFrame()
        low_expression_group = pd.DataFrame()
        for tissue in df['tissue'].unique():
            tissue_df = df[df['tissue'] == tissue]
            median = tissue_df[genes].quantile(0.5)
            high_expression_group = pd.concat([high_expression_group, tissue_df[(tissue_df[genes] > median).T.all()]])
            low_expression_group = pd.concat([low_expression_group, tissue_df[(tissue_df[genes] <= median).T.all()]])

    # Count number of patients
    n_patients_low = low_expression_group.shape[0]
    n_patients_high = high_expression_group.shape[0]
    
    try:
        df['group'] = None
        df.loc[low_expression_group.index, 'group'] = 0
        df.loc[high_expression_group.index, 'group'] = 1

        cols = ['OS.time', 'OS', 'group']
        if pan_cancer:
            cols.append('tissue')
        df = df[cols]
        df = df.dropna()
        
        # Fit Cox Proportional Hazards model
        cph = CoxPHFitter(penalizer=0.1, l1_ratio=0.)
        strata = 'tissue' if pan_cancer else None
        cph = cph.fit(df, duration_col='OS.time', event_col='OS', strata=strata)
        
        # Get Hazard Ratio
        hr = cph.hazard_ratios_['group']

        # Get p-value
        pval = cph.summary.loc['group', 'p']

        # Get AIC
        AIC = cph.AIC_partial_

        # Get concordance index
        c = cph.concordance_index_

        # p-value of the interaction, run logrank to compare the two groups
        #logrank_pval = logrank_test(
        #    durations_A=df[df['group'] == 0]['OS.time'],
        #    durations_B=df[df['group'] == 1]['OS.time'],
        #    event_observed_A=df[df['group'] == 0]['OS'],
        #    event_observed_B=df[df['group'] == 1]['OS']
        #).p_value

        # Return metrics and data for plotting
        return {
            #'pval': logrank_pval,
            'pval': pval,
            'hr': hr,
            'AIC': AIC,
            'concordance': c,
            'n_patients_low': n_patients_low,
            'n_patients_high': n_patients_high,
        }

    except Exception as e:
        print(f'Error: {e}')
        print('Interaction:', interactions)
        print('n_patients_low:', n_patients_low)
        print('n_patients_high:', n_patients_high)

        return {
            'pval': np.nan,
            'hr': np.nan,
            'AIC': np.nan,
            'concordance': np.nan,
            'n_patients_low': n_patients_low,
            'n_patients_high': n_patients_high,
        }

# Test all ccc interactions
ccc_results = pd.DataFrame(index=ccc, columns=['interaction', 'pval', 'hazard_ratio', 'AIC', 'concordance', 'n_patients_low', 'n_patients_high'])
for interaction in ccc:
    result = survival(interaction, df)
    ccc_results.loc[interaction] = { 
        'interaction': interaction,
        'pval': result['pval'],
        'hazard_ratio': result['hr'],
        'AIC': result['AIC'],
        'concordance': result['concordance'],
        'n_patients_low': result['n_patients_low'],
        'n_patients_high': result['n_patients_high'],
    }

ccc_results['type'] = 'ccc'

# Also add crosstalk interactions from the motifs file
crosstalk_results = pd.DataFrame(index=motifs, columns=['interaction', 'pval', 'hazard_ratio', 'AIC', 'concordance', 'n_patients_low', 'n_patients_high'])
for interaction in motifs:
    result = survival(interaction, df)
    crosstalk_results.loc[interaction] = { 
        'interaction': interaction,
        'pval': result['pval'],
        'hazard_ratio': result['hr'],
        'AIC': result['AIC'],
        'concordance': result['concordance'],
        'n_patients_low': result['n_patients_low'],
        'n_patients_high': result['n_patients_high'],
    }

crosstalk_results['type'] = 'crosstalk'

# Combine the results
all_results = pd.concat([ccc_results, crosstalk_results])
all_results = all_results.sort_values('pval')

# Save results
output_dir = '/home/lnemati/pathway_crosstalk/results/survival/tissues'
tissue_name = os.path.basename(args.tissuefile).split('.')[0]
output_file = os.path.join(output_dir, f'{tissue_name}.csv')

print('Saving results to:', output_file)
all_results.to_csv(output_file, index=False)

print('Done: survival_all.py')
