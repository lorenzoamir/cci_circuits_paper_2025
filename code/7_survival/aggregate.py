import os
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--inputdir', type=str, help='Path to folder contining tissue-specific results')
parser.add_argument('--outputdir', type=str, help='Path to folder to save aggregated results')

args = parser.parse_args()

path = args.inputdir
dfs = []
rawdfs = []

MIN_PATIENTS = 10

def extract_values(col):  
    df[col]     = rawdf.loc[pairs, col].values
    df[col+'1'] = rawdf.loc[int1,       col].values
    df[col+'2'] = rawdf.loc[int2,       col].values

for file in os.listdir(path):
    rawdf = pd.read_csv(os.path.join(path, file), index_col='interaction')
    #rawdf = rawdf.dropna()
    tissue = file.replace('.csv', '').title()
    print(tissue)

    rawdf['tissue'] = tissue
    rawdfs.append(rawdf)
    
    pairs = rawdf.query('(type == "crosstalk")').index
    
    if len(pairs) == 0:
        continue
    
    int1 = pairs.str.split('&', expand=True).get_level_values(0)
    int2 = pairs.str.split('&', expand=True).get_level_values(1)
    
    df = pd.DataFrame(index=pairs)

    cols = [
        'tissue',
        'hr',
        'n_patients_low',
        'n_patients_high',
        'concordance_index',
        'logrank_pval',
        'se',
        #'ci_low',
        #'ci_high',
    ]
    
    for col in cols:
        extract_values(col)
    
    df = df.drop(columns=['tissue1', 'tissue2', 'se1', 'se2'])
    
    df['min_patients']  = df[['n_patients_low', 'n_patients_high']].min(1)
    df['min_patients1'] = df[['n_patients_low1', 'n_patients_high2']].min(1)
    df['min_patients2'] = df[['n_patients_low1', 'n_patients_high2']].min(1)
    
    df = df.drop(columns=['n_patients_low', 'n_patients_low1', 'n_patients_low2', 'n_patients_high', 'n_patients_high1', 'n_patients_high2'])
    
    dfs.append(df)

rawdf =  pd.concat(rawdfs).reset_index().set_index(['interaction', 'tissue'])
pairdf = pd.concat(dfs).reset_index().set_index(['interaction', 'tissue'])

pairdf['log2_hr'] = np.log2(pairdf['hr'])
pairdf['log2_hr1'] = np.log2(pairdf['hr1'])
pairdf['log2_hr2'] = np.log2(pairdf['hr2'])

print('Filtering min patients')
valid = pairdf.dropna().query('(min_patients > @MIN_PATIENTS) and (min_patients1 > @MIN_PATIENTS) and (min_patients2 > @MIN_PATIENTS)').index

print('Filtering p-values')
rawdf.loc[valid, 'pval_adj'] = false_discovery_control(rawdf.loc[valid, 'logrank_pval'])
pairdf.loc[valid, 'pval_adj'] = rawdf.loc[valid, 'pval_adj']
pairdf = pairdf.loc[valid]

pairdf = pairdf.sort_values(by='pval_adj')
pairdf['hr_lfc_best'] =  pairdf['log2_hr'].abs() - pairdf[['log2_hr1', 'log2_hr2']].abs().max(1)
pairdf['c_diff_best'] =  pairdf['concordance_index'].abs() - pairdf[['concordance_index1', 'concordance_index2']].abs().max(1)

rawdf = rawdf.reset_index()
pairdf = pairdf.reset_index()
rawdf.set_index('interaction', inplace=True)
pairdf.set_index('interaction', inplace=True)

motifs = pd.read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv', index_col='Interaction')
rawdf.loc[rawdf.type == 'ccc', 'motif'] = 'ccc'
rawdf.loc[rawdf.index.intersection(motifs.index), 'motif'] = motifs.loc[rawdf.index.intersection(motifs.index), 'Type']
pairdf['motif'] = motifs.loc[pairdf.index, 'Type']

print('Saving')
outdir = args.outputdir
os.makedirs(outdir, exist_ok=True)
rawdf.to_csv(os.path.join(outdir, 'all_unfiltered.csv'))
pairdf.to_csv(os.path.join(outdir, 'all_pairs.csv'))

print('Finding better pairs')
better = pairdf.copy()
better = better[better.concordance_index >= 0.5]
better = better.query('pval_adj < 0.05')
better = better.query('(logrank_pval < logrank_pval1) and (logrank_pval < logrank_pval2)')
better = better.query('hr_lfc_best > 0.')
better = better.query('c_diff_best > 0.')

better = better.sort_values(by='c_diff_best', ascending=False)
better.to_csv(os.path.join(outdir, 'better_pairs.csv'))
print('Done: aggregate.py')
